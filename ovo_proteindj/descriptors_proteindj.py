from ovo.core.database.models import Descriptor, NumericGlobalDescriptor, StructureFileDescriptor, FileDescriptor

# Files

# fold_id', 'seq_id', 'rfd_sampled_mask', 'rfd_helices', 'rfd_strands', 'rfd_total_ss', 'rfd_RoG',
# 'mpnn_score', 'af2_pae_interaction', 'af2_pae_overall', 'af2_pae_binder', 'af2_pae_target', 'af2_plddt_overall',
# 'af2_plddt_binder', 'af2_plddt_target', 'af2_rmsd_overall', 'af2_rmsd_binder_bndaln', 'af2_rmsd_binder_tgtaln',
# 'af2_rmsd_target', 'pr_helices', 'pr_strands', 'pr_total_ss', 'pr_RoG', 'pr_intface_BSA', 'pr_intface_shpcomp',
# 'pr_intface_hbonds', 'pr_intface_deltaG', 'pr_intface_packstat', 'pr_TEM', 'pr_surfhphobics_%', 'seq_ext_coef',
# 'seq_length', 'seq_MW', 'seq_pI', 'sequence', 'rfd_time', 'af2_time

RFDIFFUSION_STRUCTURE_PATH = StructureFileDescriptor(
    name="RFdiffusion backbone design",
    description="RFdiffusion-generated backbone structure, without any side-chains, with Glycine residues at designed positions",
    tool="RFdiffusion",
    key="proteindj|rfd|backbone_structure_path",
    structure_type="backbone_design",
)

RFDIFFUSION_TRB_PATH = FileDescriptor(
    name="RFdiffusion .trb file",
    description="RFdiffusion pickle file containing metadata about the designed backbone structure",
    tool="RFdiffusion",
    key="proteindj|rfd|backbone_trb_path",
)

PROTEINDJ_PROTEIN_MPNN_STRUCTURE_PATH = StructureFileDescriptor(
    name="ProteinMPNN sequence design",
    description="Structure for ProteinMPNN-generated sequence",
    tool="ProteinMPNN",
    key="proteindj|mpnn|proteinmpnn_structure_path",
    structure_type="sequence_design",
)

PROTEINDJ_AF2_STRUCTURE_PATH = StructureFileDescriptor(
    name="AlphaFold2 structure prediction",
    description="AlphaFold2-predicted structure for designed sequence",
    tool="AlphaFold2",
    key="proteindj|af2|alphafold_structure_path",
    structure_type="prediction",
)

FOLD_HELICES = NumericGlobalDescriptor(
    name="Number of helices",
    description="Number of alpha-helices in the RFdiffusion-generated fold.",
    tool="RFdiffusion",
    key="proteindj|fold|fold_helices",
)

FOLD_STRANDS = NumericGlobalDescriptor(
    name="Number of strands",
    description="Number of beta-strands in the RFdiffusion-generated fold.",
    tool="RFdiffusion",
    key="proteindj|fold|fold_strands",
)

FOLD_TOTAL_SS = NumericGlobalDescriptor(
    name="Number of secondary structures",
    description="Total secondary structures (helices + strands).",
    tool="RFdiffusion",
    key="proteindj|fold|fold_total_ss",
)

FOLD_ROG = NumericGlobalDescriptor(
    name="RFdiffusion radius of gyration",
    description="Radius of gyration for RFdiffusion-generated fold.",
    tool="RFdiffusion",
    key="proteindj|fold|fold_RoG",
)

RFD_TIME = NumericGlobalDescriptor(
    name="RFdiffusion time",
    description="Computation time in seconds for RFdiffusion.",
    tool="RFdiffusion",
    key="proteindj|rfd|rfd_time",
    comparison="lower_is_better",
)

FAMPNN_AVG_PSCE = NumericGlobalDescriptor(
    name="FAMPNN average PSCE",
    description="Average predicted sidechain error; lower is better.",
    tool="FAMPNN",
    key="proteindj|mpnn|fampnn_avg_psce",
    comparison="lower_is_better",
)

MPNN_SCORE = NumericGlobalDescriptor(
    name="ProteinMPNN score",
    description="Average negative log likelihood from ProteinMPNN; lower is better.",
    tool="ProteinMPNN",
    key="proteindj|mpnn|mpnn_score",
    comparison="lower_is_better",
)

AF2_PAE_INTERACTION = NumericGlobalDescriptor(
    name="AF2 interaction PAE",
    description="Predicted aligned error for interaction interfaces; lower is better.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_pae_interaction",
    color_scale="pae",
    comparison="lower_is_better",
    min_value=0,
)

AF2_PAE_OVERALL = NumericGlobalDescriptor(
    name="AF2 overall PAE",
    description="Predicted aligned error across all chains.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_pae_overall",
    color_scale="pae",
    comparison="lower_is_better",
    min_value=0,
)

AF2_PAE_BINDER = NumericGlobalDescriptor(
    name="AF2 binder PAE",
    description="Predicted aligned error for the binder chain.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_pae_binder",
    color_scale="pae",
    comparison="lower_is_better",
    min_value=0,
)

AF2_PAE_TARGET = NumericGlobalDescriptor(
    name="AF2 target PAE",
    description="Predicted aligned error for the target.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_pae_target",
    color_scale="pae",
    comparison="lower_is_better",
    min_value=0,
)

AF2_PLDDT_OVERALL = NumericGlobalDescriptor(
    name="AF2 overall pLDDT",
    description="Average per-residue pLDDT (0–100).",
    tool="AlphaFold2",
    key="proteindj|af2|af2_plddt_overall",
    color_scale="plddt",
    min_value=0,
    max_value=100,
    comparison="higher_is_better",
)

AF2_PLDDT_BINDER = NumericGlobalDescriptor(
    name="AF2 binder pLDDT",
    description="Binder pLDDT confidence.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_plddt_binder",
    comparison="higher_is_better",
    color_scale="plddt",
    min_value=0,
    max_value=100,
)

AF2_PLDDT_TARGET = NumericGlobalDescriptor(
    name="AF2 target pLDDT",
    description="Target pLDDT confidence.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_plddt_target",
    color_scale="plddt",
    min_value=0,
    max_value=100,
    comparison="higher_is_better",
)

AF2_RMSD_OVERALL = NumericGlobalDescriptor(
    name="AF2 overall RMSD",
    description="Cα RMSD between design and AF2 prediction (all chains).",
    tool="AlphaFold2",
    key="proteindj|af2|af2_rmsd_overall",
    comparison="lower_is_better",
    color_scale="rmsd",
)

AF2_RMSD_BINDER_BNDALN = NumericGlobalDescriptor(
    name="AF2 Binder-aligned Binder RMSD",
    description="Binder RMSD between AF2 prediction and design when binder chains are aligned.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_rmsd_binder_bndaln",
    comparison="lower_is_better",
    color_scale="rmsd",
)

AF2_RMSD_BINDER_TGTALN = NumericGlobalDescriptor(
    name="AF2 Target-aligned Binder RMSD",
    description="Binder RMSD between AF2 prediction and design when target chains are aligned.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_rmsd_binder_tgtaln",
    comparison="lower_is_better",
    color_scale="rmsd",
)

AF2_RMSD_TARGET = NumericGlobalDescriptor(
    name="AF2 target RMSD",
    description="Target-chain Cα RMSD between design and AF2 model.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_rmsd_target",
    comparison="lower_is_better",
    color_scale="rmsd",
)

AF2_TIME = NumericGlobalDescriptor(
    name="AF2 inference time",
    description="Time (seconds) for AF2 prediction.",
    tool="AlphaFold2",
    key="proteindj|af2|af2_time",
    comparison="lower_is_better",
)

BOLTZ_OVERALL_RMSD = NumericGlobalDescriptor(
    name="Boltz RMSD overall",
    description="Cα RMSD of Boltz-2 prediction vs input design.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_overall_rmsd",
    comparison="lower_is_better",
    color_scale="rmsd",
)

BOLTZ_BINDER_RMSD = NumericGlobalDescriptor(
    name="Boltz RMSD binder",
    description="Binder Cα RMSD vs input design.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_binder_rmsd",
    comparison="lower_is_better",
    color_scale="rmsd",
)

BOLTZ_TARGET_RMSD = NumericGlobalDescriptor(
    name="Boltz RMSD target",
    description="Target Cα RMSD vs input design.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_target_rmsd",
    comparison="lower_is_better",
    color_scale="rmsd",
)

BOLTZ_CONF_SCORE = NumericGlobalDescriptor(
    name="Boltz confidence score",
    description="Model confidence; lower is more confident.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_confidence_score",
    comparison="lower_is_better",
    color_scale="rmsd",
)

BOLTZ_PTM = NumericGlobalDescriptor(
    name="Boltz pTM",
    description="Predicted TM-score (0–1).",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_ptm",
    comparison="higher_is_better",
)

BOLTZ_PTM_INTERFACE = NumericGlobalDescriptor(
    name="Boltz pTM interface",
    description="Interface TM-score (0–1).",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_ptm_interface",
    comparison="higher_is_better",
)

BOLTZ_PLDDT = NumericGlobalDescriptor(
    name="Boltz pLDDT",
    description="Predicted LDDT score for complex (0–1).",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_plddt",
    comparison="higher_is_better",
    color_scale="plddt",
)

BOLTZ_PLDDT_INTERFACE = NumericGlobalDescriptor(
    name="Boltz interface pLDDT",
    description="Interface pLDDT score.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_plddt_interface",
    comparison="higher_is_better",
    color_scale="plddt",
)

BOLTZ_PDE = NumericGlobalDescriptor(
    name="Boltz PDE",
    description="Predicted distance error; lower is better.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_pde",
    comparison="lower_is_better",
    color_scale="pae",
)

BOLTZ_PDE_INTERFACE = NumericGlobalDescriptor(
    name="Boltz PDE interface",
    description="Predicted distance error at interface.",
    tool="Boltz-2",
    key="proteindj|boltz|boltz_pde_interface",
    comparison="lower_is_better",
    color_scale="pae",
)

PR_HELICES = NumericGlobalDescriptor(
    name="Num helices in pred. structure",
    description="Number of alpha-helices in predicted structure.",
    tool="PyRosetta",
    key="proteindj|pr|pr_helices",
    comparison="higher_is_better",
)

PR_STRANDS = NumericGlobalDescriptor(
    name="Num strands in pred. structure",
    description="Number of beta-strands in predicted structure.",
    tool="PyRosetta",
    key="proteindj|pr|pr_strands",
    comparison="higher_is_better",
)

PR_TOTAL_SS = NumericGlobalDescriptor(
    name="Num secondary structures in pred. structure",
    description="Total helices + strands.",
    tool="PyRosetta",
    key="proteindj|pr|pr_total_ss",
    comparison="higher_is_better",
)

PR_ROG = NumericGlobalDescriptor(
    name="Radius of gyration of pred. structure",
    description="RoG for predicted structure.",
    tool="PyRosetta",
    key="proteindj|pr|pr_RoG",
    comparison="lower_is_better",
)

PR_INTERFACE_BSA = NumericGlobalDescriptor(
    name="Interface buried surface area",
    description="Buried surface area at interface (Å²).",
    tool="PyRosetta",
    key="proteindj|pr|pr_intface_BSA",
    comparison="higher_is_better",
)

PR_INTERFACE_SHPCOMP = NumericGlobalDescriptor(
    name="Interface shape complementarity",
    description="Shape complementarity (0–1).",
    tool="PyRosetta",
    key="proteindj|pr|pr_intface_shpcomp",
    comparison="higher_is_better",
)

PR_INTERFACE_HBONDS = NumericGlobalDescriptor(
    name="Interface hydrogen bonds",
    description="Number of interface hydrogen bonds.",
    tool="PyRosetta",
    key="proteindj|pr|pr_intface_hbonds",
    comparison="higher_is_better",
)

PR_INTERFACE_DELTAG = NumericGlobalDescriptor(
    name="Interface dG",
    description="Solvation free energy gain (kcal/mol).",
    tool="PyRosetta",
    key="proteindj|pr|pr_intface_deltaG",
    comparison="lower_is_better",
)

PR_INTERFACE_PACKSTAT = NumericGlobalDescriptor(
    name="Interface packstat",
    description="Interface packing quality (0–1).",
    tool="PyRosetta",
    key="proteindj|pr|pr_intface_packstat",
    comparison="higher_is_better",
)

PR_TEM = NumericGlobalDescriptor(
    name="Total energy metric",
    description="Rosetta-style total energy; lower is more stable.",
    tool="PyRosetta",
    key="proteindj|pr|pr_TEM",
    comparison="lower_is_better",
)

PR_SURF_HPHOBIC_PCT = NumericGlobalDescriptor(
    name="Surface hydrophobicity %",
    description="Percentage of hydrophobic residues exposed.",
    tool="PyRosetta",
    key="proteindj|pr|pr_surfhphobics_%",
    comparison="lower_is_better",
)

SEQ_EXT_COEF = NumericGlobalDescriptor(
    name="Extinction coefficient",
    description="Extinction coefficient at 280 nm.",
    tool="SequenceAnalysis",
    key="proteindj|seq|seq_ext_coef",
)

SEQ_LENGTH = NumericGlobalDescriptor(
    name="Sequence length",
    description="Number of residues.",
    tool="SequenceAnalysis",
    key="proteindj|seq|seq_length",
)

SEQ_MW = NumericGlobalDescriptor(
    name="Molecular weight",
    description="Sequence molecular weight in Daltons.",
    tool="SequenceAnalysis",
    key="proteindj|seq|seq_MW",
)

SEQ_PI = NumericGlobalDescriptor(
    name="Isoelectric point",
    description="Predicted isoelectric point.",
    tool="SequenceAnalysis",
    key="proteindj|seq|seq_pI",
)

DESCRIPTORS = [v for v in globals().values() if isinstance(v, Descriptor)]

PRESETS = [
    {
        "label": "AF2 iPAE & RMSD",
        "x": AF2_PAE_INTERACTION,
        "y": AF2_RMSD_BINDER_TGTALN,
    },
    {
        "label": "AF2 PAE & RMSD",
        "x": AF2_PAE_OVERALL,
        "y": AF2_RMSD_OVERALL,
    },
    {
        "label": "Boltz iPDE & RMSD",
        "x": BOLTZ_PDE_INTERFACE,
        "y": BOLTZ_BINDER_RMSD,
    },
]
