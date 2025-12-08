import streamlit as st

def proteindj_logo():
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
