import streamlit as st

from ovo import (
    local_scheduler,
    storage,
    config,
)
from ovo.app.components import molstar_custom_component, StructureVisualization
from ovo.app.components.history_components import history_dropdown_component
from ovo.app.components.input_components import pdb_input_component, sequence_selection_fragment
from ovo.app.components.molstar_custom_component import ContigsParser
from ovo.app.components.navigation import show_prev_next_sections
from ovo.app.components.preview_components import visualize_rfdiffusion_preview, contigs_organizer_fragment
from ovo.app.components.scheduler_components import wait_with_statusbar
from ovo.app.components.submission_components import (
    pool_submission_inputs,
    review_workflow_submission,
)
from ovo.app.pages.rfdiffusion.binder_diversification import initialize_workflow
from ovo.app.utils.page_init import initialize_page
from ovo.core.logic.design_logic_rfdiffusion import submit_rfdiffusion_preview
from ovo.core.utils.formatting import get_hashed_path_for_bytes
from ovo_proteindj.models_proteindj import ProteinDJMonomerMotifScaffDesignWorkflow


@st.fragment
def intro_step():
    # Initialize the workflow object in session state
    initialize_workflow(__file__, ProteinDJMonomerMotifScaffDesignWorkflow.name)

    if not config.props.pyrosetta_license:
        return st.error("This workflow requires a PyRosetta license which is not enabled in this OVO instance. Please contact the administrator.")

    with st.container(width=850):
        st.markdown(
            f"""

```
██████╗ ██████╗  ██████╗ ████████╗███████╗██╗███╗   ██╗██████╗      ██╗
██╔══██╗██╔══██╗██╔═══██╗╚══██╔══╝██╔════╝██║████╗  ██║██╔══██╗     ██║
██████╔╝██████╔╝██║   ██║   ██║   █████╗  ██║██╔██╗ ██║██║  ██║     ██║
██╔═══╝ ██╔══██╗██║   ██║   ██║   ██╔══╝  ██║██║╚██╗██║██║  ██║██   ██║
██║     ██║  ██║╚██████╔╝   ██║   ███████╗██║██║ ╚████║██████╔╝╚█████╔╝
╚═╝     ╚═╝  ╚═╝ ╚═════╝    ╚═╝   ╚══════╝╚═╝╚═╝  ╚═══╝╚═════╝  ╚════╝ 
                   ProteinDJ Protein Design Pipeline                   
          Developers: Dylan Silke, Josh Hardy, Julie Iskander       
               https://github.com/PapenfussLab/proteindj   
```

This workflow enables designing a scaffold to support segments from an existing input structure. 

TODO readme.
""",
            unsafe_allow_html=True,
        )


@st.fragment
def input_step():
    workflow: ProteinDJMonomerMotifScaffDesignWorkflow = st.session_state.workflows[__file__]

    left, _, right = st.columns([4, 1, 3])
    with left:
        st.subheader("Input structure")
        if new_pdb_input := pdb_input_component(workflow.get_input_name()):
            input_name, pdb_input_bytes = new_pdb_input

            workflow.params.rfd_input_pdb = storage.store_file_str(
                pdb_input_bytes.decode(),
                f"project/{st.session_state.project.id}/inputs/{get_hashed_path_for_bytes(pdb_input_bytes)}/{input_name}.pdb",
                overwrite=False,
            )
            if workflow.get_contig():
                st.warning("Structure was changed, clearing contig")
                workflow.set_contig("")

    with right:
        history_dropdown_component(__file__, workflow_name=ProteinDJMonomerMotifScaffDesignWorkflow.name)

    if not workflow.get_input_pdb_path():
        # No input structure provided yet, stop here
        return

    st.subheader(workflow.get_input_name())

    molstar_custom_component(
        structures=[
            StructureVisualization(pdb=storage.read_file_str(workflow.get_input_pdb_path()), color="chain-id")
        ],
        key="input_structure",
        width=700,
        height=400,
    )


@st.fragment()
def selection_step():
    workflow: ProteinDJMonomerMotifScaffDesignWorkflow = st.session_state.workflows[__file__]

    st.subheader("Select residues to keep")

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    sequence_selection_fragment(__file__, workflow.get_input_name())


@st.fragment()
def contig_preview_step():
    workflow: ProteinDJMonomerMotifScaffDesignWorkflow = st.session_state.workflows[__file__]
    st.subheader("Enter contig and generate preview")

    st.write(
        """
        The contig defines which parts of the structure will be preserved,
        in which order they will be connected, and how many residues will be designed in between.
        """
    )

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    # Input parameters for the preview
    contigs_organizer_fragment(__file__, pdb_input_string=storage.read_file_str(workflow.get_input_pdb_path()))

    # Generate preview
    st.write("#### Generate preview")

    help_column = st.columns([2, 1])[0]

    with st.columns(3)[0]:
        if "preview_timesteps" not in st.session_state:
            st.session_state.preview_timesteps = 5
        new_timesteps = st.slider(
            "Num RFdiffusion timesteps (T)",
            min_value=1,
            max_value=20,
            value=st.session_state.preview_timesteps,
            key="timesteps_input",
        )
        if new_timesteps and new_timesteps != st.session_state.preview_timesteps:
            st.session_state.preview_timesteps = new_timesteps
            # Clear previous preview if settings changed
            workflow.preview_job_id = None

    with help_column:
        st.write(f"""
        Generate a quick RFdiffusion preview of the design with reduced number of timesteps
        ({st.session_state.preview_timesteps}/50) to verify your inputs. This step is optional.

        This should take from 30 seconds to a few minutes depending on the length of the protein.
        """)

    if st.button(":material/wand_stars: Generate preview"):
        workflow.preview_job_id = submit_rfdiffusion_preview(workflow, timesteps=st.session_state.preview_timesteps)

    # Check if needed parameters are set
    if not workflow.preview_job_id:
        return

    # Wait until task is done
    result = wait_with_statusbar(local_scheduler, workflow.preview_job_id, label="Running RFdiffusion...")

    # Exit if job failed
    if not result:
        return

    visualize_rfdiffusion_preview(workflow, output_dir=local_scheduler.get_output_dir(workflow.preview_job_id))


@st.fragment()
def inpainting_step():
    workflow: ProteinDJMonomerMotifScaffDesignWorkflow = st.session_state.workflows[__file__]

    st.subheader("Sequence inpainting")

    st.write(
        "Select input structure residues that will be kept in the structure but redesigned with ProteinMPNN. This step is optional."
    )

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    if not workflow.get_contig():
        st.error("Please provide a contig in the previous step.")
        return

    parsed_contig = ContigsParser().parse_contigs_str(workflow.get_contig())
    fixed_segments = [seg for seg in parsed_contig if seg.type == "fixed"]

    sequence_selection_fragment(__file__, workflow.get_input_name(), fixed_segments=fixed_segments, inpainting=True)


@st.fragment()
def settings_step():
    st.subheader("Settings")

    workflow: ProteinDJMonomerMotifScaffDesignWorkflow = st.session_state.workflows[__file__]

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    if not workflow.get_contig():
        st.error("Contig has not been computed. Please fill in all inputs in previous step.")
        return

    pool_submission_inputs(__file__)

    contig = st.text_input(
        "Contig",
        placeholder="A123-456/10/A567-890/...",
        value=workflow.get_contig(),
        key="input_contig",
    )

    if check_contig_parsed(contig):
        workflow.set_contig(contig)

    # TODO more advanced settings


@st.fragment
def review_step():
    st.subheader("Review settings")

    review_workflow_submission(__file__)


def check_contig_parsed(contig: str | None, verbose: bool = False) -> bool:
    parser = ContigsParser()
    if not len(contig.split(" ")) == 1:
        if verbose:
            st.warning("Support of multi-chain scaffold design is experimental, proceed with caution")

    if contig.islower():
        st.warning("Input contigs are lowercase. Workflow can behave unexpectedly.")

    try:
        _ = parser.parse_contigs_str(contig) if contig else None
        return True
    except Exception as e:
        if verbose:
            st.error(
                f"There is an error in the contig specification. Please make sure you provide valid contigs. Error: {e}"
            )
        return False

if __name__ == "__main__":
    initialize_page(page_title="ProteinDJ Scaffold design")

    show_prev_next_sections(
        key=__file__,
        title="🏗️ ProteinDJ Scaffold design",
        sections={
            "intro": intro_step,
            "input": input_step,
            "selection": selection_step,
            "contig": contig_preview_step,
            "sequence inpainting": inpainting_step,
            "settings": settings_step,
            "review": review_step,
        },
    )
