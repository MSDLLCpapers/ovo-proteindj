import os
import tempfile

import streamlit as st
from ovo.core.logic.round_logic import get_or_create_project_rounds
from ovo.app.components.submission_components import pool_submission_inputs, get_pool_inputs
from ovo.app.utils.page_init import initialize_page

from ovo_proteindj.components.help_components import proteindj_logo
from ovo_proteindj.logic import import_workflow_results

if __name__ == "__main__":
    initialize_page(page_title="ProteinDJ Import")

    st.title("🪩 ProteinDJ Import")

    proteindj_logo()
    st.markdown(
                f"""
    Import existing ProteinDJ results into OVO.
    """,
                unsafe_allow_html=True,
            )

    pool_submission_inputs(__file__)
    round_id, pool_name, pool_description = get_pool_inputs(__file__)

    uploaded_files = st.file_uploader(
        "Upload ProteinDJ result directory",
        accept_multiple_files="directory",
        key="proteindj_uploader",
    )

    design_mode = st.selectbox(
        "Design mode",
        options=["monomer_motifscaff", "binder_denovo"],
        index=None
    )

    help = None
    disabled = False
    if not uploaded_files:
        help = "No files selected"
        disabled = True
    elif not design_mode:
        help = "Please select a design mode"
        disabled = True

    if st.button("Import", type="primary", disabled=disabled, help=help):

        with tempfile.TemporaryDirectory() as tmpdirname:
            for file in uploaded_files:
                os.makedirs(os.path.join(tmpdirname, os.path.dirname(file.name)), exist_ok=True)
                with open(os.path.join(tmpdirname, file.name), "wb") as f:
                    f.write(file.getbuffer())

            dirs = os.listdir(tmpdirname)
            if not dirs:
                st.error("No files found in the uploaded directory.")
                st.stop()
            if len(dirs) > 1:
                st.error("Please upload a single ProteinDJ result directory.")
                st.stop()

            pool = import_workflow_results(
                design_mode=design_mode,
                source_dir=os.path.join(tmpdirname, dirs[0]),
                round_id=round_id,
                pool_name=pool_name,
                pool_description=pool_description,
            )
            st.success(f"Imported pool: **{pool.name}**")
