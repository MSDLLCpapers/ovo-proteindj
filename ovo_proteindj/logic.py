import glob
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from typing import Callable

import pandas as pd
from ovo import (
    db,
    storage,
    get_scheduler,
    Design,
    get_username,
    Scheduler,
)
from ovo.core.database import descriptors

from ovo.core.database.models import Pool, Round, DesignJob, DesignSpec, DescriptorValue, Base, WorkflowTypes
from ovo.core.logic.descriptor_logic import save_descriptor_job_for_design_job, read_descriptor_file_values, \
    generate_descriptor_values_for_design
from ovo.core.logic.design_logic import set_designs_accepted

from ovo_proteindj import models_proteindj
from ovo_proteindj import descriptors_proteindj


def process_workflow_results(job: DesignJob, callback: Callable = None) -> list[Base]:
    pool = db.get(Pool, design_job_id=job.id)
    scheduler = get_scheduler(job.scheduler_key)
    assert job.workflow.is_instance(models_proteindj.ProteinDJDesignWorkflow), (
        "This function expects a ProteinDJDesignWorkflow instance, got: {}".format(type(job.workflow).__name__)
    )
    return process_proteindj_output_dir(
        job=job,
        pool=pool,
        source_dir=scheduler.get_output_dir(job.job_id),
        fold_dir="rfd",
        seq_dir=job.workflow.params.seq_method,
        pred_dir=job.workflow.params.pred_method,
        callback=callback,
    )


def import_workflow_results(design_mode: str, source_dir: str, round_id: str, pool_name: str, pool_description: str = ""):
    """Import ProteinDJ workflow results from a local directory into a given project round.

    Saves Pool, DesignJob, Designs and DescriptorValues into the database.
    """
    source_dir = os.path.abspath(source_dir)
    assert os.path.exists(source_dir), "Source directory does not exist: {}".format(source_dir)
    assert os.path.isdir(source_dir), "Source directory is not a directory: {}".format(source_dir)
    if db.count(Pool, round_id=round_id, name=pool_name):
        raise ValueError(f"Pool with name '{pool_name}' already exists in round {round_id}")
    # Create pool and dummy job
    project_round = db.get(Round, id=round_id)
    pool = Pool(
        id=Pool.generate_id(),
        round_id=round_id,
        name=pool_name,
        description=pool_description,
        author=get_username()
    )
    WorkflowSubclass = WorkflowTypes.get(f"ProteinDJ {design_mode}")
    input_filenames = sorted(glob.glob(os.path.join(source_dir, "inputs/*.pdb")))
    if not input_filenames:
        raise ValueError(f"No input PDB files found in {os.path.join(source_dir, 'inputs')}")
    if len(input_filenames) > 1:
        # TODO can we have multiple input PDBs?
        raise ValueError(f"Multiple input PDB files found in {os.path.join(source_dir, 'inputs')}, expected only one, got: {input_filenames}")

    input_filename = os.path.basename(input_filenames[0])

    job = DesignJob(
        scheduler_key=Scheduler.IMPORTED_SCHEDULER_KEY,
        job_id=Scheduler.IMPORTED_JOB_ID,
        workflow=WorkflowSubclass(
            params={
                "input_pdb": storage.store_input(
                    project_id=project_round.project_id,
                    file_path=os.path.join(source_dir, "inputs", input_filename),
                ),
                # TODO read params from source_dir
            },
            # TODO read thresholds from source_dir, or get them from the user
            #  or modify processing logic to ensure that acceptance is based on membership in best_designs.csv
            # acceptance_thresholds=...
        ),
        job_result=True,
        author=get_username(),
    )
    # Detect backbone, sequence design and structure prediction directories
    run_dir = os.path.join(source_dir, "run")
    fold_dir = None
    for option in ["rfd"]:
        if os.path.exists(os.path.join(source_dir, f"run/{option}/")):
            fold_dir = option
            break
    assert fold_dir is not None, (f"Could not find backbone generation directory in {run_dir}, "
                                     f"found run dirs: {', '.join(os.listdir(run_dir))}")
    seq_dir = None
    for option in ["fampnn", "mpnn"]:
        if os.path.exists(os.path.join(source_dir, f"run/{option}/")):
            seq_dir = option
            break
    assert seq_dir is not None, (f"Could not find sequence design directory in {run_dir}, "
                                 f"found run dirs: {', '.join(os.listdir(run_dir))}")
    pred_dir = None
    for option in ["af2", "boltz"]:
        if os.path.exists(os.path.join(source_dir, f"run/{option}/")):
            pred_dir = option
            break
    assert pred_dir is not None, (f"Could not find structure prediction directory in {run_dir}, "
                                 f"found run dirs: {', '.join(os.listdir(run_dir))}")
    # Process the output directory
    objects = [pool, job] + process_proteindj_output_dir(
        job=job,
        pool=pool,
        source_dir=source_dir,
        fold_dir=fold_dir,
        seq_dir=seq_dir,
        pred_dir=pred_dir,
    )
    db.save_all(objects)
    # Update design_job_id in pool
    pool.design_job_id = job.id
    db.save(pool)
    return pool


def process_proteindj_output_dir(
    job: DesignJob,
    pool: Pool,
    source_dir: str,
    fold_dir: str,
    seq_dir: str,
    pred_dir: str,
    callback: Callable = None
) -> list[Base]:

    # this is where result files will be stored in our storage
    # make sure to remove trailing slash otherwise S3 will keep two slashes in the path
    project_round = db.get(Round, id=pool.round_id)
    destination_dir = os.path.join("project", project_round.project_id, "pools", pool.id, "designs").rstrip("/")

    all_designs_df = pd.read_csv(StringIO(storage.read_file_str(os.path.join(source_dir, "results/all_designs.csv"))), index_col=0)
    assert all_designs_df.index.is_unique, "Design IDs in all_designs.csv are not unique: {}".format(all_designs_df.index.tolist())

    num_fold_designs = all_designs_df.fold_id.max() + 1
    num_sequence_designs = all_designs_df.seq_id.max() + 1

    # TODO we should use a temp directory near our storage to speed up file operations below
    with tempfile.TemporaryDirectory() as tmpdir:
        for dirname in [fold_dir, seq_dir, pred_dir]:
            # sync the result file (when using S3 storage, this downloads the file, otherwise we just symlink it)
            storage.sync_file(
                os.path.join(source_dir, f"run/{dirname}/{dirname}_results.tar.gz"),
                os.path.join(tmpdir, f"{dirname}_results.tar.gz"),
                link=True
            )
            # untar the result files
            tar_path = os.path.join(tmpdir, f'{dirname}_results.tar.gz')
            subprocess.run(['tar', '-xzf', tar_path, '-C', tmpdir], check=True)

        designs = []
        descriptor_values = []
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(
                    process_proteindj_design,
                    pool_id=pool.id,
                    destination_dir=destination_dir,
                    tmpdir=tmpdir,
                    fold_dir=fold_dir,
                    seq_dir=seq_dir,
                    pred_dir=pred_dir,
                    num_fold_designs=num_fold_designs,
                    num_sequence_designs=num_sequence_designs,
                    row=row
                )
                for _, row in all_designs_df.iterrows()
            ]

            for i, future in enumerate(futures):
                new_design, new_descriptors = future.result()
                designs.append(new_design)
                descriptor_values.extend(new_descriptors)
                if callback:
                    callback(
                        value=(i + 1) / len(futures),
                        text=f"Downloading design {new_design.id}",
                    )

    # Update design.accepted fields based on descriptor values and thresholds
    set_designs_accepted(designs, descriptor_values, job.workflow.acceptance_thresholds)

    # Save descriptor job and assign its ID to descriptor values
    descriptor_job = save_descriptor_job_for_design_job(
        design_job=job,
        project_id=project_round.project_id,
        chains=["FULL"],
        design_ids=[design.id for design in designs],
    )
    for descriptor_value in descriptor_values:
        descriptor_value.descriptor_job_id = descriptor_job.id

    return designs + descriptor_values


def process_proteindj_design(
    pool_id: str,
    destination_dir: str,
    tmpdir: str,
    fold_dir: str,
    seq_dir: str,
    pred_dir: str,
    num_fold_designs: int,
    num_sequence_designs: int,
    row: pd.Series
) -> tuple[Design, list[DescriptorValue]]:
    total_idx_fold = row.fold_id
    idx_sequence = row.seq_id
    fold_filename = f"fold_{total_idx_fold}"

    # fold suffix 01, 001, 0001 based on total number of designs
    fold_suffix = "_" + str(total_idx_fold + 1).zfill(max(len(str(num_fold_designs)), 2))
    # 01, 001 based on total number of sequences
    sequence_id = str(idx_sequence + 1).zfill(max(len(str(num_sequence_designs)), 2))

    fold_id = f"ovo_{pool_id}{fold_suffix}"
    design_id = fold_id + f"_seq{sequence_id}"

    # create Design object
    design = Design(
        id=design_id,
        pool_id=pool_id,
        accepted=False,
    )
    rfdiffusion_backbone_pdb_path = storage.store_file_path(
        source_abs_path=f"{tmpdir}/{fold_dir}_results/{fold_filename}.pdb",
        storage_rel_path=f"{destination_dir}/rfdiffusion/{fold_id}_backbone.pdb",
        overwrite=False,
    )

    fold_json_path = storage.store_file_path(
        source_abs_path=os.path.join(tmpdir, f"{fold_dir}_results/{fold_filename}.json"),
        storage_rel_path=f"{destination_dir}/rfdiffusion/{fold_id}.json",
        overwrite=False,
    )
    sequence_design_descriptor = descriptors_proteindj.PROTEINDJ_PROTEIN_MPNN_STRUCTURE_PATH
    sequence_design_pdb_path = storage.store_file_path(
        source_abs_path=f"{tmpdir}/{seq_dir}_results/{fold_filename}_seq_{idx_sequence}.pdb",
        storage_rel_path=f"{destination_dir}/mpnn/{design_id}.pdb",
        overwrite=False,
    )
    #
    # TODO handle Boltz prediction here and in generate_descriptor_values_for_design below
    #
    structure_prediction_descriptor = descriptors_proteindj.PROTEINDJ_AF2_STRUCTURE_PATH
    structure_prediction_pdb_path = storage.store_file_path(
        source_abs_path=f"{tmpdir}/{pred_dir}_results/{fold_filename}_seq_{idx_sequence}_af2pred.pdb",
        storage_rel_path=f"{destination_dir}/alphafold_initial_guess/{design_id}_af2pred.pdb",
        overwrite=False,
    )
    shared_args = dict(
        design_id=design.id,
        descriptor_job_id=None,
        chains="A",
    )
    descriptor_values = [
        DescriptorValue(
            descriptor_key=descriptors_proteindj.RFDIFFUSION_STRUCTURE_PATH.key,
            value=rfdiffusion_backbone_pdb_path,
            **shared_args,
        ),
        DescriptorValue(
            descriptor_key=descriptors_proteindj.FOLD_JSON_PATH.key,
            value=fold_json_path,
            **shared_args,
        ),
        DescriptorValue(
            descriptor_key=sequence_design_descriptor.key,
            value=sequence_design_pdb_path,
            **shared_args,
        ),
        DescriptorValue(
            descriptor_key=structure_prediction_descriptor.key,
            value=structure_prediction_pdb_path,
            **shared_args,
        ),
    ]
    # Use the prediction PDB as the design structure - it's guaranteed to have valid sidechains
    design.structure_path = structure_prediction_pdb_path
    # TODO handle multiple chains
    chain_ids = ["A"]
    design.spec = DesignSpec.from_pdb_str(
        pdb_data=storage.read_file_str(design.structure_path), chains=chain_ids
    )
    # TODO populate each spec.chains[...].contig from RFD json or csv row,
    #  make sure they are in the correct order as the chains in the PDB!

    # Sanity check - make sure sequences in CSV match those in PDB
    if ":" in row.sequence:
        spec_sequences = "|".join([f"{chain_id}:{c.sequence}" for c in design.spec.chains for chain_id in c.chain_ids])
    else:
        assert len(design.spec.chains), f"Expected one chain in design {design.id}, got {len(design.spec.chains)}"
        spec_sequences = design.spec.chains[0].sequence
    assert spec_sequences == str(row.sequence), \
        f"Unexpected mismatch between sequence in PDB and CSV for design {design.id}: {spec_sequences} != {row.sequence}"

    df = pd.DataFrame([row.rename(design_id)])

    # Legacy support for rfd_ prefix instead of fold_
    df = df.rename(columns={
        "rfd_helices": "fold_helices",
        "rfd_strands": "fold_strands",
        "rfd_total_ss": "fold_total_ss",
        "rfd_RoG": "fold_RoG",
    })

    descriptor_values.extend(
        generate_descriptor_values_for_design(
            design_id=design_id,
            table_ids=design_id,
            descriptor_job_id=None, # will be set later
            descriptor_tables={
                "proteindj|fold": df,
                "proteindj|rfd": df,
                "proteindj|mpnn": df,
                "proteindj|af2": df,
                # TODO include boltz metrics
                #"proteindj|boltz": df,
                "proteindj|pr": df,
                "proteindj|seq": df,
            },
            chains=chain_ids,
        )
    )

    return design, descriptor_values
