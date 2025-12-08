import streamlit as st
from ovo_proteindj.models_proteindj import default_params


def proteindj_intro():
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
    """,
        unsafe_allow_html=True,
    )

    if not default_params:
        st.warning(
            "ProteinDJ default parameters are not configured. Please add the following block to your OVO config file:"
        )
        st.code(
            """
plugins:
 ovo_proteindj:
   default_params:
     rfd_models: "{config.reference_files_dir}/rfdiffusion_models"
     af2_models: "{config.reference_files_dir}/alphafold_models"
     boltz_models: "{config.reference_files_dir}/boltz_models"
     cpus: 4
     cpus_per_gpu: 4
     memory_gpu: 14GB
     memory_cpu: 14GB
     gpu_queue: gpu-small
     gpu_model: a10g
            """,
            language="yaml",
        )
