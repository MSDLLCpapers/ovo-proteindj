import datetime
import os

import pytest

from ovo.core.utils.resources import RESOURCES_DIR
from ovo.core.utils.tests import TEST_PROJECT_NAME


@pytest.fixture(scope="session", autouse=True)
def check_unit_test_mode():
    if os.environ.get("OVO_UNIT_TEST") == "1":
        raise RuntimeError("Please run unit tests and workflow tests separately.")


@pytest.fixture(scope="session")
def project_data():
    """Create one project and project_round for the entire test run."""
    from ovo import db, storage, Design, Pool

    # Use repository root directory name as round prefix
    prefix = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    round_name = prefix + " " + datetime.datetime.now().strftime("%a %d %b %Y")

    if not db.Project.count(name=TEST_PROJECT_NAME):
        project = db.Project(name=TEST_PROJECT_NAME, author="test", public=True)
        db.save(project)
    else:
        project = db.Project.get(name=TEST_PROJECT_NAME)

    if db.Round.count(name=round_name, project_id=project.id):
        # Clear round if exists
        project_round = db.Round.get(name=round_name, project_id=project.id)
        pools = db.Pool.select(round_id=project_round.id)
        designs = db.Design.select(pool_id__in=[p.id for p in pools])
        db.DescriptorValue.remove(design_id__in=[d.id for d in designs])
        db.Design.remove(pool_id__in=[p.id for p in pools])
        db.DesignJob.remove(id__in=[p.design_job_id for p in pools if p.design_job_id])
        db.Pool.remove(round_id=project_round.id)
        db.Round.remove(id=project_round.id)

    project_round = db.Round(
        name=round_name,
        project_id=project.id,
        author="test",
    )
    db.save(project_round)

    # Create a pool
    custom_pool = Pool(
        id=Pool.generate_id(),
        round_id=project_round.id,
        name="Custom Upload",
        author="test",
    )
    db.save(custom_pool)

    # Create a test design
    test_design = Design.from_pdb_file(
        storage=storage,
        filename="5ELI_A.pdb",
        pdb_str=(RESOURCES_DIR / "examples/inputs/5ELI_A.pdb").read_text(),
        chains=["A"],
        project_id=project.id,
        pool_id=custom_pool.id,
    )
    db.save(test_design)
    assert len(test_design.spec.chains) == 1
    assert os.path.exists(os.path.join(storage.storage_root, test_design.structure_path))

    return project, project_round, custom_pool
