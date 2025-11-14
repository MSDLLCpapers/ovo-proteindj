import json
import os
import pickle
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
)
from ovo.core.database import descriptors

from ovo.core.database.models import Pool, Round, DesignJob, DesignSpec, DescriptorValue
from ovo.core.logic.descriptor_logic import save_descriptor_job_for_design_job
from ovo.core.logic.design_logic import set_designs_accepted
from ovo.core.scheduler.base_scheduler import Scheduler

from ovo_proteindj.models import ProteinDJDesignWorkflow


def submit_workflow(workflow: ProteinDJDesignWorkflow, scheduler: Scheduler, pipeline_name: str = None):
    workflow.validate()

    # TODO workflow input files should be prepared automatically by the scheduler in the future
    params = workflow.params.to_dict()
    params["rfd_input_pdb"] = storage.prepare_workflow_input(params["rfd_input_pdb"], scheduler.workdir)

    job_id = scheduler.submit(
        pipeline_name=pipeline_name or "https://github.com/PapenfussLab/proteindj@v1.0.0",
        params=params,
        # TODO do we need to pass -profile binder_denovo to set default values,
        #  or will we set all default values on our end?
    )

    return job_id


def process_workflow_results(job: DesignJob, callback: Callable = None):
    pool = db.get(Pool, design_job_id=job.id)
    project_round = db.get(Round, id=pool.round_id)
    scheduler = get_scheduler(job.scheduler_key)
    workflow: ProteinDJDesignWorkflow = job.workflow
    assert workflow.is_instance(ProteinDJDesignWorkflow), (
        "This function expects a ProteinDJDesignWorkflow instance, got: {}".format(type(workflow).__name__)
    )

    # this is where result files will be stored in our storage
    # make sure to remove trailing slash otherwise S3 will keep two slashes in the path
    destination_dir = os.path.join("project", project_round.project_id, "pools", pool.id, "designs").rstrip("/")

    source_output_path = scheduler.get_output_dir(job.job_id)

    all_designs_df = pd.read_csv(StringIO(storage.read_file_str(os.path.join(source_output_path, "results/all_designs.csv"))), index_col=0)
    assert all_designs_df.index.is_unique, "Design IDs in all_designs.csv are not unique: {}".format(all_designs_df.index.tolist())

    num_backbone_designs = all_designs_df.fold_id.max() + 1
    num_sequence_designs = all_designs_df.seq_id.max() + 1

    rfd_dir = "rfd"
    seq_dir = workflow.params.seq_method
    pred_dir = workflow.params.pred_method

    with tempfile.TemporaryDirectory() as tmpdir:
        for dirname in [rfd_dir, seq_dir, pred_dir]:
            # untar the result files
            storage.sync_file(os.path.join(source_output_path, f"run/{dirname}/{dirname}_results.tar.gz"), tmpdir)
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
                    rfd_dir=rfd_dir,
                    seq_dir=seq_dir,
                    pred_dir=pred_dir,
                    num_backbone_designs=num_backbone_designs,
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
    set_designs_accepted(designs, descriptor_values, workflow.acceptance_thresholds)
    # Save descriptor job and assign its ID to descriptor values
    descriptor_job = save_descriptor_job_for_design_job(
        design_job=job,
        round_id=pool.round_id,
        chains=["A"],
        design_ids=[design.id for design in designs],
    )
    for descriptor_value in descriptor_values:
        descriptor_value.descriptor_job_id = descriptor_job.id
    # Save designs and descriptors atomically in one commit
    db.save_all(designs + descriptor_values)


def process_proteindj_design(
    pool_id: str,
    destination_dir: str,
    tmpdir: str,
    rfd_dir: str,
    seq_dir: str,
    pred_dir: str,
    num_backbone_designs: int,
    num_sequence_designs: int,
    row: pd.Series
) -> tuple[Design, list[DescriptorValue]]:
    total_idx_backbone = row.fold_id
    idx_sequence = row.seq_id
    backbone_filename = f"fold_{total_idx_backbone}"

    # backbone suffix 01, 001, 0001 based on total number of designs
    backbone_suffix = "_" + str(total_idx_backbone + 1).zfill(max(len(str(num_backbone_designs)), 2))
    # 01, 001 based on total number of sequences
    sequence_id = str(idx_sequence + 1).zfill(max(len(str(num_sequence_designs)), 2))

    backbone_id = f"ovo_{pool_id}{backbone_suffix}"
    design_id = backbone_id + f"_seq{sequence_id}"

    rfd_json = json.loads(storage.read_file_str(
        os.path.join(tmpdir, f"{rfd_dir}_results/{backbone_filename}.json")
    ))

    # create Design object
    design = Design(
        id=design_id,
        pool_id=pool_id,
        accepted=False,
        rfdiffusion_backbone_pdb_path=storage.store_file_path(
            source_abs_path=f"{tmpdir}/{rfd_dir}_results/{backbone_filename}.pdb",
            storage_rel_path=f"{destination_dir}/rfdiffusion/{backbone_id}_backbone.pdb",
            overwrite=False,
        ),
        rfdiffusion_backbone_trb_path=storage.store_file_bytes(
            file_bytes=pickle.dumps(rfd_json),
            storage_rel_path=f"{destination_dir}/rfdiffusion/{backbone_id}_backbone.trb",
            overwrite=False,
        ),
        protein_mpnn_pdb_path=storage.store_file_path(
            source_abs_path=f"{tmpdir}/{seq_dir}_results/{backbone_filename}_seq_{idx_sequence}.pdb",
            storage_rel_path=f"{destination_dir}/mpnn/{design_id}.pdb",
            overwrite=False,
        ),
        alphafold2_initial_guess_pdb_path=storage.store_file_path(
            source_abs_path=f"{tmpdir}/{pred_dir}_results/{backbone_filename}_seq_{idx_sequence}_af2pred.pdb",
            storage_rel_path=f"{destination_dir}/alphafold_initial_guess/{design_id}_af2pred.pdb",
            overwrite=False,
        )
    )
    design.structure_path = design.protein_mpnn_pdb_path
    chain_ids = ["A"]
    design.spec = DesignSpec.from_pdb_str(
        pdb_data=storage.read_file_str(design.structure_path), chains=chain_ids
    )

    # Sanity check - make sure sequences in CSV match those in PDB
    if ":" in row.sequence:
        spec_sequences = "|".join([f"{chain_id}:{c.sequence}" for c in design.spec.chains for chain_id in c.chain_ids])
    else:
        assert len(design.spec.chains), f"Expected one chain in design {design.id}, got {len(design.spec.chains)}"
        spec_sequences = design.spec.chains[0].sequence
    assert spec_sequences == str(row.sequence), \
        f"Unexpected mismatch between sequence in PDB and CSV for design {design.id}: {spec_sequences} != {row.sequence}"

    #
    # TODO should we share descriptors or rather create new ones?
    #
    mapping = {
        # Unique Identifiers
        "description": None,
        "fold_id": None,
        "seq_id": None,

        # RFdiffusion Metrics
        "rfd_sampled_mask": None,
        "rfd_helices": None,
        "rfd_strands": None,
        "rfd_total_ss": None,
        "rfd_RoG": descriptors.RADIUS_OF_GYRATION.key,
        "rfd_time": None,

        # Sequence Generation Metrics
        "fampnn_avg_psce": None,
        "mpnn_score": None,

        # Structure Prediction Metrics
        "af2_pae_interaction": descriptors.AF2_BINDER_INTERFACE_PAE.key, # TODO rename to binder_interaction_pae
        "af2_pae_overall": descriptors.AF2_PAE.key,
        "af2_pae_binder": descriptors.AF2_BINDER_PAE.key,
        "af2_pae_target": None,
        "af2_plddt_overall": descriptors.AF2_PLDDT.key,
        "af2_plddt_binder": descriptors.AF2_PLDDT_BINDER.key,
        "af2_plddt_target": None,
        "af2_rmsd_overall": descriptors.AF2_DESIGN_RMSD.key,
        "af2_rmsd_binder_bndaln": None,
        "af2_rmsd_binder_tgtaln": descriptors.AF2_TARGET_ALIGNED_BINDER_RMSD.key,
        "af2_rmsd_target": None,
        "af2_time": None,

        "boltz_overall_rmsd": None, # TODO add boltz metrics
        "boltz_binder_rmsd": None,
        "boltz_target_rmsd": None,
        "boltz_conf_score": None,
        "boltz_ptm": None,
        "boltz_ptm_interface": None,
        "boltz_plddt": None,
        "boltz_plddt_interface": None,
        "boltz_pde": None,
        "boltz_pde_interface": None,

        # Biophysical Analysis Metrics
        "pr_helices": None,
        "pr_strands": None,
        "pr_total_ss": None,
        "pr_RoG": None,
        "pr_intface_BSA": None,
        "pr_intface_shpcomp": None,
        "pr_intface_hbonds": None,
        "pr_intface_deltaG": descriptors.PYROSETTA_DDG.key,
        "pr_intface_packstat": None,
        "pr_TEM": None,
        "pr_surfhphobics_%": None,

        # Sequence Analysis Metrics
        "seq_ext_coef": None,
        "seq_length": descriptors.LENGTH.key,
        "seq_MW": descriptors.MOLECULEAR_WEIGHT.key,
        "seq_pI": descriptors.ISOELECTRIC_POINT.key,
    }

    descriptor_values = []
    for column, value in row.items():
        if column in mapping and mapping[column]:
            assert isinstance(mapping[column], str), f"Expected descriptor key, got {type(mapping[column])}: {mapping[column]}"
            descriptor_values.append(
                DescriptorValue(
                    design_id=design.id,
                    descriptor_key=mapping[column],
                    descriptor_job_id=None, # will be set later
                    chains=",".join(chain_ids),
                    value=str(value) if pd.notna(value) else None,
                )
            )

    return design, descriptor_values
