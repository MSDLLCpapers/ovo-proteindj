from ovo import db
from ovo.core.logic import design_logic
from ovo.core.utils.tests import TEST_SCHEDULER_KEY
from ovo_proteindj import models_proteindj, descriptors_proteindj
from ovo.core.utils.resources import RESOURCES_DIR

def test_binder_denovo(project_data):
    project, project_round, custom_pool = project_data

    workflow = models_proteindj.ProteinDJBinderDeNovoDesignWorkflow(
        params=models_proteindj.ProteinDJBinderDeNovoDesignWorkflow.Params(
            input_pdb=RESOURCES_DIR / "examples/inputs/5ELI_A.pdb",
            rfd_contigs="[A74-97/0 20]",
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
        pool_name="5ELI ProteinDJ binder_denovo",
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

    af2_ipae = db.select_descriptor_values(descriptors_proteindj.AF2_PAE_INTERACTION.key, design_ids)
    assert len(af2_ipae.dropna()) == 2
    assert (af2_ipae < 30).all()

    rosetta_ddg = db.select_descriptor_values(descriptors_proteindj.PR_INTERFACE_DELTAG.key, design_ids)
    assert len(rosetta_ddg.dropna()) == 2
    assert not rosetta_ddg.isna().any()
