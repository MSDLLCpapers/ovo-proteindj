import json
import os

import streamlit as st

from ovo import db, storage
from ovo.app.components.molstar_custom_component import (
    molstar_custom_component,
    StructureVisualization,
    ChainVisualization,
    ContigsParser,
)
from ovo.app.utils.cached_db import (
    get_cached_design_descriptors,
    get_cached_design,
    get_cached_pool,
    get_cached_design_job,
    get_cached_designs,
)
from ovo.core.database import (
    Design,
    descriptors,
)
from ovo.core.database.descriptors import (
    SEQUENCE_DESIGN_PATH_DESCRIPTORS,
)
from ovo.core.database.models_rfdiffusion import (
    RFdiffusionWorkflow,
    RFdiffusionBinderDesignWorkflow,
    RFdiffusionScaffoldDesignWorkflow,
)
from ovo.core.utils.colors import get_color_from_str
from ovo.core.utils.pdb import align_multiple_proteins_pdb, pdb_to_mmcif, filter_pdb_str, get_sequences_from_pdb_str
from ovo.core.utils.residue_selection import from_segments_to_hotspots, from_hotspots_to_segments
from ovo.app.components.workflow_visualization_components import visualize_design_sequence, show_design_metrics, visualize_scaffold_alignment

from ovo_proteindj import descriptors_proteindj
import re


def parse_pdb_idx(idx_str: str) -> list[int]:
    if isinstance(idx_str, list):
        return idx_str
    assert isinstance(idx_str, str), f"Index must be a string or list of strings, got {type(idx_str)}: {idx_str}"
    # TODO remove this logic once the idx is stored as proper json
    # replace np.int64(56) -> 56
    idx_str = re.sub(r"np\.int\d+\((\d+)\)", r"\1", idx_str)
    idx_str = idx_str.replace("'", '"').replace("(", "[").replace(")", "]")
    return json.loads(idx_str)


def get_output_segments(parser, fold_meta: dict):
    return parser.parse_contigs_ref(
        contigs=fold_meta["rfd_sampled_mask"][0],
        ref_idx=parse_pdb_idx(fold_meta["rfd_con_ref_pdb_idx"]),
        hal_idx=parse_pdb_idx(fold_meta["rfd_con_hal_pdb_idx"]),
    )


def visualize_scaffold_design_structure(design_id: str | None):
    design: Design = get_cached_design(design_id)

    # TODO show contigs once they are in spec
    #chain_contigs = [c.contig + "/0" for c in design.spec.chains]
    #st.write(f"Contig: **{' '.join(chain_contigs)}**")

    show_design_metrics(
        design_id,
        descriptor_keys=[
            descriptors_proteindj.SEQ_LENGTH.key,
            descriptors_proteindj.FOLD_ROG.key,
            descriptors_proteindj.FOLD_HELICES.key,
            descriptors_proteindj.FOLD_STRANDS.key,
            descriptors_proteindj.AF2_PAE_OVERALL.key,
            descriptors_proteindj.AF2_RMSD_OVERALL.key,
            # TODO include motif rmsd when available
            #descriptors_proteindj.AF2_RMSD_MOTIF.key,
            descriptors_proteindj.AF2_PLDDT_OVERALL.key,
        ],
    )

    pool = get_cached_pool(design.pool_id)
    design_job = get_cached_design_job(pool.design_job_id)
    workflow: RFdiffusionScaffoldDesignWorkflow = design_job.workflow

    parser = ContigsParser()
    paths = (
        get_cached_design_descriptors(
            design_id,
            [
                descriptors_proteindj.FOLD_JSON_PATH.key,
                *[d.key for d in descriptors.STRUCTURE_PATH_DESCRIPTORS],
            ],
        )
        .dropna()
        .to_dict()
    )
    fold_meta = json.loads(storage.read_file_str(paths[descriptors_proteindj.FOLD_JSON_PATH.key]))
    input_pdb_str = storage.read_file_str(workflow.get_input_pdb_path())
    input_segments = [segment for contig in fold_meta["rfd_sampled_mask"] for segment in parser.parse_contigs_str(contig)]

    output_segments = get_output_segments(parser, fold_meta)
    input_mapping = [
        (segment.input_res_chain, list(range(segment.input_res_start, segment.input_res_end + 1)))
        for segment in input_segments
        if segment.type == "fixed"
    ]
    output_mapping = [
        (segment.out_res_chain, list(range(segment.out_res_start, segment.out_res_end + 1)))
        for segment in output_segments
        if segment.type == "fixed"
    ]

    left, middle, right = st.columns(3, gap="medium")

    with left:
        st.write("##### Full input structure")

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=input_pdb_str, contigs=input_segments, representation_type="cartoon+ball-and-stick"
                ),
            ],
            key="full_input",
            height=350,
        )

    backbone_design_descriptor = descriptors_proteindj.RFDIFFUSION_STRUCTURE_PATH

    with middle:
        st.write(f"##### {backbone_design_descriptor.name}")

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=storage.read_file_str(paths[backbone_design_descriptor.key]),
                    contigs=output_segments,
                )
            ],
            key="rfdiff_contig_segments",
            height=350,
        )

        st.write(backbone_design_descriptor.description)

    sequence_design_descriptor = None
    for d in SEQUENCE_DESIGN_PATH_DESCRIPTORS:
        if paths.get(d.key):
            sequence_design_descriptor = d

    with right:
        st.write(f"##### {sequence_design_descriptor.name}")

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=storage.read_file_str(paths[sequence_design_descriptor.key]),
                    contigs=output_segments,
                    representation_type="cartoon+ball-and-stick",
                )
            ],
            key="mpnn_contig_segments",
            height=350,
        )

        st.write(sequence_design_descriptor.description)

    st.write("### Refolding tests")

    prediction_descriptor = descriptors_proteindj.PROTEINDJ_AF2_STRUCTURE_PATH

    if not prediction_descriptor:
        return

    left, middle, right = st.columns(3, gap="medium")

    prediction_storage_path = paths[prediction_descriptor.key]
    prediction_pdb = storage.read_file_str(prediction_storage_path)

    with left:
        st.write(f"##### {prediction_descriptor.name}")
        description = prediction_descriptor.description

        if prediction_descriptor.b_factor_value == "plddt":
            pdb_str = pdb_to_mmcif(prediction_pdb, "-", True)
            description += " (colored by pLDDT confidence)"
        elif prediction_descriptor.b_factor_value == "fractional_plddt":
            pdb_str = pdb_to_mmcif(prediction_pdb, "-", True, fractional_plddt=True)
            description += " (colored by pLDDT confidence)"
        else:
            pdb_str = prediction_pdb

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=pdb_str,
                    representation_type="cartoon",
                    color="plddt",
                ),
            ],
            key=f"pred_{prediction_descriptor.key}",
            height=350,
        )
        st.write(description)

    with middle:
        st.write(f"##### Input motif aligned to prediction")

        input_motif_pdb = filter_pdb_str(input_pdb_str, [s.value for s in input_segments if s.type == "fixed"])
        structures = [
            StructureVisualization(
                pdb=input_motif_pdb, contigs=input_segments, representation_type="cartoon+ball-and-stick"
            )
        ]

        (_, aligned_prediction_pdb), rmsd = align_multiple_proteins_pdb(
            pdb_strs=[input_motif_pdb, prediction_pdb],
            chain_residue_mappings=[
                input_mapping,
                output_mapping,
            ],
            all_atom=True,
        )
        num_fixed = sum(
            segment.out_res_end - segment.out_res_start + 1 for segment in output_segments if segment.type == "fixed"
        )

        structures.append(
            StructureVisualization(
                pdb=aligned_prediction_pdb,
                representation_type="cartoon+ball-and-stick",
            )
        )

        molstar_custom_component(
            structures=structures,
            key=f"input_motif_aligned_to_{prediction_descriptor.name}",
            height=350,
            html_filename=os.path.basename(prediction_storage_path).replace(".pdb", "_aligned_motif"),
        )

        st.write(
            f"""
            {prediction_descriptor.name} vs native motif all atom RMSD: **{rmsd:.2f} Å**

            {prediction_descriptor.name} prediction (white) aligned to original input structure (colored segments) 
            of the {num_fixed} fixed input motif residues.
            """
        )

    with right:
        st.write(f"##### Design aligned to prediction")

        structures, rmsd = align_multiple_proteins_pdb(
            pdb_strs=[
                storage.read_file_str(paths[sequence_design_descriptor.key]),
                prediction_pdb,
            ],
            chain_residue_mappings=[[("A", None)], [("A", None)]],
            all_atom=False,
        )
        aligned_str = structures[1]

        if prediction_descriptor.b_factor_value == "plddt":
            aligned_str = pdb_to_mmcif(aligned_str, "-", True)
        elif prediction_descriptor.b_factor_value == "fractional_plddt":
            aligned_str = pdb_to_mmcif(aligned_str, "-", True, fractional_plddt=True)

        molstar_custom_component(
            structures=[
                StructureVisualization(pdb=structures[0], contigs=output_segments),
                StructureVisualization(
                    pdb=aligned_str,
                    representation_type="cartoon",
                    color="plddt",
                ),
            ],
            key=f"design_aligned_to_{prediction_descriptor.key}",
            height=350,
        )

        st.write(
            f"""
            {prediction_descriptor.name} vs design backbone RMSD: **{rmsd:.2f} Å**

            Agreement between {sequence_design_descriptor.name} and {prediction_descriptor.name}
            """
        )


def visualize_scaffold_design_sequence(design_id: str):
    # visualize regular sequence
    visualize_design_sequence(design_id)

    # visualize alignment
    st.write("##### Alignment of input sequence to designed sequence")

    design: Design = get_cached_design(design_id)
    pool = get_cached_pool(design.pool_id)
    workflow: RFdiffusionWorkflow = get_cached_design_job(pool.design_job_id).workflow

    input_pdb_path = workflow.get_input_pdb_path()
    input_seq_by_resno: dict[str, dict[str, str]] = get_sequences_from_pdb_str(
        storage.read_file_str(input_pdb_path), by_residue_number=True
    )
    designed_sequences: dict[str, str] = {
        chain_id: chain.sequence for chain in design.spec.chains for chain_id in chain.chain_ids
    }
    paths = (
        get_cached_design_descriptors(design_id, [descriptors_proteindj.FOLD_JSON_PATH.key]).dropna().to_dict()
    )

    parser = ContigsParser()
    # TODO this could be parsed from spec contig instead
    fold_meta = json.loads(
        storage.read_file_str(paths[descriptors_proteindj.FOLD_JSON_PATH.key])
    )
    output_segments = get_output_segments(parser, fold_meta)
    # individual positions in this format ['A1', 'A2', ...]
    inpainted_positions = (
        from_segments_to_hotspots(workflow.get_inpaint_seq().split("/"))
        if workflow.get_inpaint_seq()
        else []
    )
    visualize_scaffold_alignment(
        input_seq_by_resno=input_seq_by_resno,
        designed_sequences=designed_sequences,
        parsed_segments=output_segments,
        inpainted_positions=inpainted_positions,
    )



def visualize_binder_design_structure(design_id: str):
    design: Design = get_cached_design(design_id)
    pool = get_cached_pool(design.pool_id)
    design_job = get_cached_design_job(pool.design_job_id)
    workflow: RFdiffusionBinderDesignWorkflow = design_job.workflow
    paths = (
        get_cached_design_descriptors(
            design_id,
            [
                *[d.key for d in descriptors.STRUCTURE_PATH_DESCRIPTORS],
            ],
        )
        .dropna()
        .to_dict()
    )
    input_pdb_path = workflow.get_input_pdb_path()

    st.write(f"Input structure: **{os.path.basename(input_pdb_path)}**")

    show_design_metrics(
        design_id,
        descriptor_keys=[
            descriptors_proteindj.SEQ_LENGTH.key,
            descriptors_proteindj.FOLD_ROG.key,
            descriptors_proteindj.FOLD_HELICES.key,
            descriptors_proteindj.FOLD_STRANDS.key,
            descriptors_proteindj.PR_INTERFACE_DELTAG.key,
            descriptors_proteindj.PR_INTERFACE_HBONDS.key,
            descriptors_proteindj.PR_INTERFACE_PACKSTAT.key,
            descriptors_proteindj.PR_INTERFACE_BSA.key,
            descriptors_proteindj.AF2_PAE_INTERACTION.key,
            descriptors_proteindj.AF2_PLDDT_BINDER.key,
            descriptors_proteindj.AF2_RMSD_BINDER_TGTALN.key,
            descriptors_proteindj.AF2_RMSD_BINDER_BNDALN.key,
        ],
    )

    # TODO use target spec for this
    # Note this assumes that the target chain is B in the RFdiffusion output PDB,
    # and that the residue numbers are same as in the input PDB
    hotspot_selections = (
        [h.replace(workflow.target_chain, "B") for h in workflow.get_hotspots().split(",")]
        if workflow.get_hotspots() and workflow.target_chain
        else None
    )

    left, middle, right = st.columns(3, gap="medium")

    backbone_design_descriptor = descriptors_proteindj.RFDIFFUSION_STRUCTURE_PATH
    with left:
        st.write(f"##### {backbone_design_descriptor.name}")

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=storage.read_file_str(paths[backbone_design_descriptor.key]),
                    highlighted_selections=hotspot_selections,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            color_params={"value": get_color_from_str(design.id, "seaborn:tab10_light")},
                            representation_type="cartoon",
                            label=f"{backbone_design_descriptor.name} {design.id}",
                        ),
                        ChainVisualization(chain_id="B", representation_type="cartoon", label=f"{design.id} target"),
                    ],
                )
            ],
            key="backbone_design",
            height=350,
        )

        st.write(backbone_design_descriptor.description)

    sequence_design_descriptor = None
    for d in SEQUENCE_DESIGN_PATH_DESCRIPTORS:
        if paths.get(d.key):
            sequence_design_descriptor = d

    with middle:
        st.write(f"##### {sequence_design_descriptor.name}")

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=storage.read_file_str(paths[sequence_design_descriptor.key]),
                    highlighted_selections=hotspot_selections,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            color_params={"value": get_color_from_str(design.id, "seaborn:tab10_light")},
                            representation_type="cartoon+ball-and-stick",
                            label=f"{sequence_design_descriptor.name} {design.id}",
                        ),
                        ChainVisualization(
                            chain_id="B", representation_type="cartoon+ball-and-stick", label=f"{design.id} target"
                        ),
                    ],
                )
            ],
            key="sequence_design_1",
            height=350,
        )

        st.write(sequence_design_descriptor.description)

    with right:
        st.write("##### Designed binding pose")

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=storage.read_file_str(paths[sequence_design_descriptor.key]),
                    highlighted_selections=hotspot_selections,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            color_params={"value": get_color_from_str(design.id, "seaborn:tab10_light")},
                            representation_type="cartoon+ball-and-stick",
                            label=f"{sequence_design_descriptor.name} {design.id}",
                        ),
                        ChainVisualization(
                            chain_id="B",
                            representation_type="molecular-surface",
                            color="hydrophobicity",
                            label=f"{design.id} target",
                        ),
                    ],
                )
            ],
            key="sequence_design_2",
            height=350,
        )

        st.write(
            f"{sequence_design_descriptor.name} against target surface colored by hydrophobicity scale :green-badge[**green** = hydroPHOBIC] :red-badge[**red** = hydroPHILIC]"
        )

    st.write("### Refolding tests")

    # TODO handle boltz prediction
    prediction_descriptor = descriptors_proteindj.PROTEINDJ_AF2_STRUCTURE_PATH

    if not prediction_descriptor:
        return

    left, middle, right = st.columns(3, gap="medium")
    with left:
        st.write("##### Input structure aligned to prediction")

        input_pdb_str = storage.read_file_str(input_pdb_path)
        prediction_str = storage.read_file_str(paths[prediction_descriptor.key])

        # here, manual alignment is needed
        structures, rmsd = align_multiple_proteins_pdb(
            pdb_strs=[input_pdb_str, prediction_str], chain_residue_mappings=[None, [("B", None)]]
        )
        aligned_str = structures[1]
        if prediction_descriptor.b_factor_value == "plddt":
            aligned_str = pdb_to_mmcif(aligned_str, "-", True)
        elif prediction_descriptor.b_factor_value == "fractional_plddt":
            aligned_str = pdb_to_mmcif(aligned_str, "-", True, fractional_plddt=True)

        hotspot_segments = from_hotspots_to_segments(workflow.get_hotspots()) if workflow.get_hotspots() else None
        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=structures[0],
                    representation_type="cartoon",
                    color="chain-id",
                    highlighted_selections=hotspot_segments,
                ),
                StructureVisualization(
                    pdb=aligned_str,
                    representation_type=None,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            representation_type="cartoon+ball-and-stick",
                            color="plddt",
                            label=prediction_descriptor.name,
                        ),
                        ChainVisualization(
                            chain_id="B",
                            representation_type="cartoon",
                            label="Prediction of target chain aligned to input target chain",
                        ),
                    ],
                ),
            ],
            key="full_input_aligned_to_prediction",
            height=350,
        )

        st.write(
            f"All chains from input structure aligned to {prediction_descriptor.name} of binder chain (ball and stick colored by pLDDT confidence) and target chain (white)."
        )

    with middle:
        st.write("##### Design aligned to prediction")

        # here, we do not need manual alignment, but we still do it
        structures, rmsd = align_multiple_proteins_pdb(
            pdb_strs=[
                storage.read_file_str(paths[sequence_design_descriptor.key]),
                storage.read_file_str(paths[prediction_descriptor.key]),
            ],
            chain_residue_mappings=[[("B", None)] for _ in range(2)],
        )

        aligned_str = structures[1]
        if prediction_descriptor.b_factor_value == "plddt":
            aligned_str = pdb_to_mmcif(aligned_str, "-", True)

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=structures[0],
                    representation_type=None,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            color_params={"value": get_color_from_str(design.id, "seaborn:tab10_light")},
                            representation_type="cartoon+ball-and-stick",
                            label=sequence_design_descriptor.name,
                        ),
                        ChainVisualization(chain_id="B", representation_type="cartoon", label="MPNN target"),
                    ],
                ),
                StructureVisualization(
                    pdb=aligned_str,
                    representation_type=None,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            representation_type="cartoon+ball-and-stick",
                            color="plddt",
                            label=prediction_descriptor.name,
                        ),
                        ChainVisualization(chain_id="B", representation_type="cartoon", label="Target chain"),
                    ],
                ),
            ],
            key="design_aligned_to_af2_2",
            height=350,
        )

        st.write(
            f"""
            {prediction_descriptor.name} vs design backbone RMSD: **{rmsd:.2f} Å**

            Agreement between {sequence_design_descriptor.name} and {prediction_descriptor.name}
            """
        )

    with right:
        st.write("##### Predicted binding pose")

        pdb_str = storage.read_file_str(paths[prediction_descriptor.key])
        if prediction_descriptor.b_factor_value == "plddt":
            pdb_str = pdb_to_mmcif(pdb_str, "-", True)

        molstar_custom_component(
            structures=[
                StructureVisualization(
                    pdb=pdb_str,
                    representation_type=None,
                    chains=[
                        ChainVisualization(
                            chain_id="A",
                            representation_type="ball-and-stick",
                            color="plddt",
                            label="Prediction of binder chain",
                        ),
                        ChainVisualization(
                            chain_id="B",
                            color="hydrophobicity",
                            representation_type="molecular-surface",
                            label="Prediction of target chain",
                        ),
                    ],
                ),
            ],
            key="binding_pose_predicted",
            height=350,
        )

        st.write(
            "Prediction of binder (ball and stick colored by pLDDT confidence) and target (surface colored by hydrophobicity scale) :green-badge[**green** = hydroPHOBIC] :red-badge[**red** = hydroPHILIC]"
        )
