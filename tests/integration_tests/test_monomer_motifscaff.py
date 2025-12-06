from ovo import db
from ovo.core.logic import design_logic
from ovo.core.utils.tests import TEST_SCHEDULER_KEY
from ovo_proteindj import models_proteindj, descriptors_proteindj
from ovo.core.utils.resources import RESOURCES_DIR

def test_monomer_motifscaff_hairpin(project_data):
    project, project_round, custom_pool = project_data

    workflow = models_proteindj.ProteinDJMonomerMotifScaffDesignWorkflow(
        params=models_proteindj.ProteinDJMonomerMotifScaffDesignWorkflow.Params(
            input_pdb=RESOURCES_DIR / "examples/inputs/5ELI_A.pdb",
            rfd_contigs="[A111-114/10/A117-119]",
            num_designs=1,
            seqs_per_design=2,
        )
    )
    workflow.validate()
    workflow.get_table_row()

    design_job, pool = design_logic.submit_design_workflow(
        # Your initialized workflow settings
        workflow=workflow,
        # Name your pool of designs
        pool_name="5ELI ProteinDJ monomer_motifscaff",
        pool_description="",
        # see schedulers for available scheduler keys
        scheduler_key=TEST_SCHEDULER_KEY,
        # Project and Round where the Pool will be created
        round_id=project_round.id,
    )

    jobs = design_logic.get_design_jobs_table(id=pool.id)
    print(jobs)
    assert len(jobs) == 1

    # wait for job to complete and process results
    pool = design_logic.process_results(design_job)

    num_designs = db.Design.count(pool_id=pool.id)
    assert num_designs == 2

    designs = db.Design.select(pool_id=pool.id)
    design_ids = [d.id for d in designs]

    num_helices = db.select_descriptor_values(descriptors_proteindj.FOLD_HELICES.key, design_ids)
    assert len(num_helices.dropna()) == 2
    assert (num_helices >= 0).all()

    af2_pae = db.select_descriptor_values(descriptors_proteindj.AF2_PAE_OVERALL.key, design_ids)
    assert len(af2_pae.dropna()) == 2
    assert (af2_pae < 30).all()

    design_rmsd = db.select_descriptor_values(descriptors_proteindj.AF2_RMSD_OVERALL, design_ids)
    assert len(design_rmsd.dropna()) == 2
    assert (design_rmsd < 20).all()

