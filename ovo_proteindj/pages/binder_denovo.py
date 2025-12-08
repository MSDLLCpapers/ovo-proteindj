import re

import streamlit as st
from ovo import config, local_scheduler, storage
from ovo.app.components.history_components import history_dropdown_component
from ovo.app.components.input_components import pdb_input_component, sequence_selection_fragment
from ovo.app.components.molstar_custom_component import molstar_custom_component, StructureVisualization
from ovo.app.components.navigation import show_prev_next_sections
from ovo.app.components.preview_components import visualize_rfdiffusion_preview
from ovo.app.components.scheduler_components import wait_with_statusbar
from ovo.app.components.submission_components import (
    pool_submission_inputs,
    review_workflow_submission,
)
from ovo.app.components.trim_components import check_hotspots, trimmed_structure_visualizer
from ovo.app.pages.rfdiffusion.binder_diversification import initialize_workflow
from ovo.app.utils.page_init import initialize_page
from ovo.core.logic.design_logic_rfdiffusion import submit_rfdiffusion_preview
from ovo.core.utils.formatting import get_hashed_path_for_bytes
from ovo.core.utils.residue_selection import from_contig_to_residues, from_residues_to_chain_breaks, \
    get_chains_and_contigs, from_residues_to_segments

from ovo_proteindj.models_proteindj import ProteinDJBinderDeNovoDesignWorkflow


@st.fragment
def intro_step():
    # Initialize the workflow object in session state
    initialize_workflow(__file__, ProteinDJBinderDeNovoDesignWorkflow.name)

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

This workflow enables designing new peptide and miniprotein binders to your input structure. 

TODO readme.
""",
            unsafe_allow_html=True,
        )


@st.fragment
def input_step():
    workflow: ProteinDJBinderDeNovoDesignWorkflow = st.session_state.workflows[__file__]

    left, _, right = st.columns([4, 1, 3])
    with left:
        st.subheader("Input structure")
        if new_pdb_input := pdb_input_component(workflow.get_input_name()):
            input_name, pdb_input_bytes = new_pdb_input

            workflow.set_input_pdb_path(
                storage.store_input(
                    project_id=st.session_state.project.id,
                    filename=f"{input_name}.pdb",
                    file_bytes=pdb_input_bytes,
                )
            )
            if workflow.get_hotspots() or workflow.get_target_contig():
                st.warning("Structure was changed, clearing hotspots and trimming region")
                workflow.set_contig("")
                workflow.set_hotspots(None)

    with right:
        history_dropdown_component(__file__, workflow_name=ProteinDJBinderDeNovoDesignWorkflow.name)

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


@st.fragment
def hotspots_step():
    workflow: ProteinDJBinderDeNovoDesignWorkflow = st.session_state.workflows[__file__]

    st.subheader("Select binding hotspots")

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    sequence_selection_fragment(__file__, workflow.get_input_name(), color="hydrophobicity", write_segments=False)



def parameters_trim_structure_component(workflow, pdb_input_string: str):
    st.write("Trim the target chain to save computation time, and to avoid running out of GPU memory")

    target_chain = workflow.get_target_chain()
    chain_contigs: dict[str, str] = get_chains_and_contigs(pdb_input_string)

    if workflow.get_hotspots():
        st.write(f"Binding hotspots: {workflow.get_hotspots()}")
        chains_in_hotspot = set([hotspot[0] for hotspot in workflow.get_hotspots().split(",")])
        if len(chains_in_hotspot) == 1 and target_chain is None:
            # pre-fill target chain if determined by selected hotspot
            target_chain = list(chains_in_hotspot)[0]

    if target_chain is None and len(chain_contigs) == 1:
        target_chain = list(chain_contigs)[0]

    target_chain_col, start_trimmed_col, end_trimmed_col, _ = st.columns([0.15, 0.15, 0.15, 0.45])
    with target_chain_col:
        chain_ids = sorted(chain_contigs.keys())
        try:
            selected_idx = chain_ids.index(target_chain) if target_chain else None
        except ValueError:
            # chain ID not present in our list
            selected_idx = None

        target_chain = st.selectbox(
            label="Target chain",
            placeholder=f"{len(chain_ids)} available chains",
            index=selected_idx,
            key="target_chain",
            options=chain_ids,
            help="Preview supports only one chain at the moment.",
        )

    if not target_chain:
        return

    # List residues from the selected chain
    residues = list(map(int, from_contig_to_residues(chain_contigs.get(target_chain))))
    if not workflow.get_target_contig():
        workflow.set_target_contig(chain_contigs[target_chain])

    start_res, end_res = workflow.get_target_trim_boundary()

    with start_trimmed_col:
        start_res = st.number_input(
            "Start residue",
            min_value=residues[0],
            max_value=residues[-1],
            value=start_res,
            key=f"start_trim_residue_{target_chain}",
            args=(
                target_chain,
                residues,
            ),
        )

        if start_res not in residues:
            start_res = min([res for res in residues if res > start_res], key=lambda res: abs(res - start_res))
            st.warning(
                "Start residue is not in the list of residues (between chain breaks). "
                f"Choosing {target_chain}{start_res} as start residue."
            )

    with end_trimmed_col:
        end_res = st.number_input(
            "End residue",
            min_value=residues[0],
            max_value=residues[-1],
            value=end_res,
            key=f"end_trim_residue_{target_chain}",
            args=(
                target_chain,
                residues,
            ),
        )

        if end_res not in residues:
            end_res = min([res for res in residues if res < end_res], key=lambda res: abs(res - end_res))
            st.warning(
                "End residue is not in the list of residues (between chain breaks). "
                f"Choosing {target_chain}{end_res} as end residue."
            )

    trimmed_residues = [res for res in residues if (res <= end_res) and (res >= start_res)]
    chain_break_segments = from_residues_to_chain_breaks(trimmed_residues)
    if chain_break_segments:
        st.info(
            f"Trimmed region {target_chain}{start_res}-{end_res} contains chain breaks: "
            f"{', '.join(chain_break_segments)}. "
            f"AlphaFold2 should detect these and insert a chain break in the prediction, "
            f"but please verify this in the predicted structure."
        )

    if not end_res > start_res:
        st.error("End residue must be greater than start residue.")
        return

    # Check if the hotspots are included in the selected range / chain
    if workflow.get_hotspots():
        check_hotspots(f"{target_chain}{start_res}-{end_res}", workflow.get_hotspots())

    workflow.set_target_contig(
        target_contig="/".join(from_residues_to_segments(
            chain_id=target_chain,
            residues=from_contig_to_residues(chain_contigs[target_chain]),
            start_res=start_res,
            end_res=end_res,
        ))
    )


@st.fragment()
def trim_step():
    workflow: ProteinDJBinderDeNovoDesignWorkflow = st.session_state.workflows[__file__]
    st.subheader("Target chain trimming")

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    pdb_input_string = storage.read_file_str(workflow.get_input_pdb_path())

    # Input parameters for the preview
    parameters_trim_structure_component(workflow, pdb_input_string)

    # Check if inputs are set
    if not workflow.get_target_contig():
        return

    # Visualize trimmed structure
    with st.container(height=630, border=False):
        trimmed_structure_visualizer(workflow, pdb_input_string)


def parameters_binder_preview_component(workflow):
    if not workflow.get_target_chain():
        st.error("Please provide a target chain in the previous step.")
        return

    st.write(f"Target chain: {workflow.get_target_chain()}")

    binder_length_col, hotspots_col, _ = st.columns([1, 1, 3])

    with binder_length_col:
        binder_length = st.text_input(
            "Binder length (min-max)",
            placeholder="For example 20-40",
            value=workflow.get_binder_contig(),
            key="binder_length",
        )
        if binder_length != workflow.get_binder_contig():
            # Reset the contig to force regeneration below
            workflow.set_binder_contig(binder_length)

    with hotspots_col:
        hotspots = st.text_input(
            "Hotspot residues (optional)",
            placeholder="For example A123,A124,A131",
            value=workflow.get_hotspots(),
            key="hotspots",
        )
        if hotspots:
            if not all(
                re.fullmatch("[A-Z][0-9]+", hotspot) for hotspot in hotspots.split(",")
            ):
                st.error("Invalid hotspots format, expected 'A123,A124,A131'")
                return
        if hotspots != workflow.get_hotspots():
            workflow.set_hotspots(hotspots)
            # Reset the preview_job_id to prevent the visualization of the previous preview
            workflow.preview_job_id = None

    if workflow.get_hotspots():
        check_hotspots(workflow.get_target_contig(), workflow.get_hotspots())


@st.fragment()
def preview_step():
    workflow: ProteinDJBinderDeNovoDesignWorkflow = st.session_state.workflows[__file__]
    st.subheader("Enter binder settings and preview design")

    if not workflow.get_input_pdb_path():
        st.error("Please provide an input structure in the input structure step.")
        return

    # Input parameters for the preview
    parameters_binder_preview_component(workflow)

    # Generate preview
    st.write("#### Generate preview")
    num_timesteps = 15
    with st.columns([2, 1])[0]:
        st.write(f"""
        Generate a quick RFdiffusion preview of the design with reduced number of timesteps 
        ({num_timesteps}/50) to verify your inputs. This step is optional.

        This should take 2-10 minutes depending on the length of the target and binder.
        """)

    if st.button(":material/wand_stars: Generate preview"):
        workflow.preview_job_id = submit_rfdiffusion_preview(workflow, timesteps=num_timesteps)

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
def settings_step():
    st.subheader("Settings")

    workflow: ProteinDJBinderDeNovoDesignWorkflow = st.session_state.workflows[__file__]

    st.session_state.current_pool_settings = pool_submission_inputs(__file__)

    workflow.set_contig(st.text_input(
        f"RFdiffusion contigs",
        value=workflow.get_contig(),
        key="rfd_contigs"
    ))

    workflow.params["num_designs"] = st.number_input(
        f"Number of RFdiffusion backbone designs",
        value=workflow.params.get("num_designs", 10),
        key="num_designs"
    )

    workflow.params["seqs_per_design"] = st.number_input(
        f"Number of sequences per backbone",
        value=workflow.params.get("seqs_per_design", 8),
        key="seqs_per_design"
    )


@st.fragment
def review_step():
    st.subheader("Review settings")

    review_workflow_submission(__file__)


if __name__ == '__main__':
    initialize_page(page_title="ProteinDJ Binder Design")

    if not config.props.pyrosetta_license:
        st.error("This workflow requires a PyRosetta license which is not enabled in this OVO instance. Please contact the administrator.")
        st.stop()

    show_prev_next_sections(
        key=__file__,
        title="🪩 ProteinDJ Binder Design",
        sections={
            'intro': intro_step,
            'input': input_step,
            'hotspots': hotspots_step,
            'trim': trim_step,
            'preview': preview_step,
            'settings': settings_step,
            'review': review_step,
        }
    )
