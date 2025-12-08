import os
import re
from dataclasses import dataclass, field
from typing import Callable, Literal, TypedDict, Optional

from ovo import config, storage
from ovo.core.database import DesignWorkflow, WorkflowTypes, WorkflowParams, Base, \
    DesignJob, Threshold
from ovo.core.database.models_refolding import RefoldingSupportedDesignWorkflow
from ovo.core.utils.residue_selection import from_segments_to_hotspots, \
    from_hotspots_to_segments, from_contig_to_residues

from ovo_proteindj import descriptors_proteindj

# Pipeline version
# tag or commit ID from https://github.com/PapenfussLab/proteindj
PIPELINE_VERSION = "5db7c4e6b9507ebfdcd20090ebb3eef544abe27c"

# Get default pipeline params from OVO config
plugin_config = config.plugins.get(__package__, {})
default_params = {
    key: value.format(config=config) if isinstance(value, str) else value
    for key, value in plugin_config.get("default_params", {}).items()
}

class ProteinDJParamsType(TypedDict):
    """Typed dictionary for ProteinDJ workflow parameters

    Used only for type hinting purposes - actual params are stored as a dict
    """
    input_pdb: str
    design_mode: str
    rfd_contigs: str
    seq_method: str
    hotspot_residues: str | None
    rfd_inpaint_seq: str | None
    # add other parameters as needed


@dataclass
class ProteinDJDesignWorkflow(DesignWorkflow):

    params: ProteinDJParamsType = field(default_factory=lambda: {
        **default_params,
        "out_dir": "output",
    })

    preview_job_id: str = None

    def get_pipeline_name(self) -> str:
        # TODO store pipeline version as a field in the workflow?
        #  it would be nice to know the version,
        #  but the problem is that nextflow currently only supports managing one single version of a pipeline at a time
        #  and also it might not be intuitive to users that we are not using the latest version
        #  when using the "reuse previous job" dropdown that initializes the workflow with all old fields
        return f"https://github.com/PapenfussLab/proteindj@{PIPELINE_VERSION}"

    def prepare_params(self, workdir: str) -> dict:
        params = self.params.copy()
        # TODO do this automatically in ovo design_logic - based on param schema "format" field?
        params["input_pdb"] = storage.prepare_workflow_input(params["input_pdb"], workdir)
        return params

    def process_results(self, job: "DesignJob", callback: Callable = None) -> list[Base]:
        """Process results of a successful workflow - download files from workdir, save DesignJob, Pool and Designs"""
        from ovo_proteindj.logic import process_workflow_results
        return process_workflow_results(job=job, callback=callback)

    def get_input_name(self) -> str | None:
        return os.path.basename(self.params["input_pdb"].rsplit(".", 1)[0]) if self.params.get("input_pdb") else None

    def get_input_pdb_path(self) -> str | None:
        return self.params.get("input_pdb")

    def set_input_pdb_path(self, pdb_path: str):
        self.params["input_pdb"] = pdb_path

    def get_selected_segments(self) -> list[str] | None:
        # implemented in subclasses
        raise NotImplementedError()

    def set_selected_segments(self, segments: list[str]):
        # implemented in subclasses
        raise NotImplementedError()

    def get_contig(self, contig_index: int = 0) -> list[str] | Literal[""]:
        """Get contig string without [] brackets for the given contig index (default 0)"""
        contigs = self.params["rfd_contigs"].strip("[]").split(",") if self.params.get("rfd_contigs") else []
        return contigs[contig_index] if contigs else ""

    def set_contig(self, contig: str):
        """Set contig string, adding [] brackets if not already present"""
        self.params["rfd_contigs"] = f"[{contig.strip('[]')}]" if contig else "[]"

    def get_cyclic_offset(self) -> bool:
        return False

    def get_hotspots(self) -> str | None:
        return None

    def get_refolding_design_paths(self, design_ids: list[str]) -> dict[str, str]:
        from ovo import db

        values = db.select_descriptor_values(descriptors_proteindj.PROTEINDJ_PROTEIN_MPNN_STRUCTURE_PATH.key, design_ids=design_ids)
        return values.dropna().to_dict()



@WorkflowTypes.register("ProteinDJ monomer_motifscaff")
@dataclass
class ProteinDJMonomerMotifScaffDesignWorkflow(ProteinDJDesignWorkflow, RefoldingSupportedDesignWorkflow):

    params: dict = field(default_factory=lambda: {
        **default_params,
        "out_dir": "output",
        "design_mode": "monomer_motifscaff",
        "seq_method": "fampnn",
    })

    # acceptance threshold values (descriptor key -> interval (min, max, enabled))
    acceptance_thresholds: dict[str, Threshold] = field(
        default_factory=lambda: {
            descriptors_proteindj.AF2_PAE_OVERALL.key: Threshold(max_value=5.0),
            descriptors_proteindj.AF2_RMSD_OVERALL.key: Threshold(max_value=2.0),
            # TODO implement native motif RMSD calculation
            # descriptors_proteindj.AF2_RMSD_MOTIF.key: Threshold(max_value=2.0),
            descriptors_proteindj.AF2_PLDDT_OVERALL.key: Threshold(min_value=80),
            descriptors_proteindj.FOLD_ROG.key: Threshold(enabled=False),
        }
    )

    @classmethod
    def visualize_single_design_structures(cls, design_id: str):
        """Visualize single design structures in Streamlit"""
        from ovo_proteindj.components.workflow_visualization_components import visualize_scaffold_design_structure

        visualize_scaffold_design_structure(design_id)

    @classmethod
    def visualize_single_design_sequences(self, design_id: str):
        from ovo_proteindj.components.workflow_visualization_components import visualize_scaffold_design_sequence

        visualize_scaffold_design_sequence(design_id)

    def get_inpaint_seq(self):
        return self.params.get("rfd_inpaint_seq", "").strip("[]") or None

    def get_selected_segments(self, inpainting=False) -> list[str] | None:
        if inpainting:
            return self.params["rfd_inpaint_seq"].strip("[]").split("/") if self.params.get("rfd_inpaint_seq") else []
        contig = self.get_contig()
        return [segment for subcontig in contig.split() for segment in subcontig.split("/") if segment and segment[0].isalpha()]

    def set_selected_segments(self, segments: list[str], inpainting=False) -> None:
        if inpainting:
            self.params["rfd_inpaint_seq"] = "[" + "/".join(segments) + "]" if segments else None
            return
        # TODO preserve generated segments, see ovo logic for re-generating contig
        # contig = self.get_contig()
        # generated_segments = [segment for subcontig in contig.split() for segment in subcontig.split("/") if segment and segment[0].isnumeric()]
        # if generated_segments:
        #     see ovo logic for re-generating contig
        # else:
        self.set_contig("/".join(segments))

    def get_refolding_design_type(self) -> str:
        return "scaffold"

    def get_refolding_native_pdb_path(self, contig_index: int) -> str:
        # In binder design, we compare with fixed input motif from input PDB
        return self.get_input_pdb_path()


@WorkflowTypes.register("ProteinDJ binder_denovo")
@dataclass
class ProteinDJBinderDeNovoDesignWorkflow(ProteinDJDesignWorkflow, RefoldingSupportedDesignWorkflow):

    params: dict = field(default_factory=lambda: {
        **default_params,
        "out_dir": "output",
        "design_mode": "binder_denovo",
        "hotspot_residues": None,
    })

    acceptance_thresholds: dict[str, Threshold] = field(
        default_factory=lambda: {
            descriptors_proteindj.AF2_PAE_INTERACTION.key: Threshold(max_value=10.0),
            descriptors_proteindj.AF2_RMSD_BINDER_TGTALN.key: Threshold(max_value=2.0),
            descriptors_proteindj.AF2_RMSD_BINDER_BNDALN.key: Threshold(max_value=1.0),
            descriptors_proteindj.AF2_PLDDT_BINDER.key: Threshold(min_value=80),
            descriptors_proteindj.PR_INTERFACE_DELTAG.key: Threshold(max_value=-30.0),
            descriptors_proteindj.AF2_PAE_BINDER.key: Threshold(max_value=5.0),
            descriptors_proteindj.FOLD_ROG.key: Threshold(enabled=False),
        }
    )

    @classmethod
    def visualize_single_design_structures(cls, design_id: str):
        """Visualize single design structures in Streamlit"""
        from ovo_proteindj.components.workflow_visualization_components import visualize_binder_design_structure

        visualize_binder_design_structure(design_id)

    # TODO adjust this for scaffold vs binder
    def get_selected_segments(self) -> list[str] | None:
        return from_hotspots_to_segments(self.params.get("hotspot_residues"))

    def set_selected_segments(self, segments: list[str]) -> None:
        # Convert from segments to hotspots
        residues = from_segments_to_hotspots(segments)
        self.set_hotspots(residues)

    def get_target_chain(self) -> str | None:
        if not self.get_contig():
            return None
        target_contig = self.get_target_contig()
        chains = [segment[0] for segment in target_contig.split("/") if segment and segment[0].isalpha()]
        unique_chains = sorted(set(chains))
        assert len(unique_chains) == 1, f"Expected exactly one unique chain in rfd_contigs, got: {unique_chains} in target contig {target_contig}"
        return unique_chains[0]

    def get_target_contig(self) -> str | Literal[""]:
        contig = self.get_contig()
        if not contig:
            return ""
        subcontigs = contig.split()
        assert len(subcontigs) == 2, f"Expected a binder chain and target chain in rfd_contigs, got: {contig}"
        fixed_subcontigs = [s for s in subcontigs if s and s[0].isalpha()]
        assert len(fixed_subcontigs) == 1, f"Expected exactly one fixed chain in contig: {contig}"
        return fixed_subcontigs[0]

    def get_target_trim_boundary(self) -> tuple[int | None, int | None]:
        target_contig = self.get_target_contig()
        if not target_contig:
            return None, None
        target_residues = from_contig_to_residues(target_contig)
        return target_residues[0], target_residues[-1]

    def get_binder_contig(self) -> str | Literal[""]:
        contig = self.get_contig()
        if not contig:
            return ""
        subcontigs = contig.split()
        assert len(subcontigs) == 2, f"Expected a binder chain and target chain in rfd_contigs, got: {contig}"
        designed_subcontigs = [s for s in subcontigs if s and not s[0].isalpha()]
        assert len(designed_subcontigs) == 1, f"Expected exactly one designed chain in contig: {contig}"
        return designed_subcontigs[0].removesuffix("/0")

    def set_target_contig(self, target_contig: str):
        binder_contig = self.get_binder_contig()
        if not binder_contig:
            binder_contig = "60-100/0"
        # only numbers 1-10 format
        assert re.fullmatch(r"[0-9]+(-[0-9]+)?(/0)?", binder_contig), \
            f"Expected binder contig to be in format 50-100 or 50, got: {binder_contig}"
        if not binder_contig.endswith("/0"):
            binder_contig += "/0"
        self.params["rfd_contigs"] = binder_contig + " " + target_contig

    def set_binder_contig(self, binder_contig: str):
        target_contig = self.get_target_contig()
        assert target_contig, "Target contig must be set before setting binder contig"
        if not all(re.fullmatch(r"^\d+(-\d+)?$", segment) for segment in binder_contig.split("/")):
            raise ValueError(f"Expected binder contig to be in format 50-100 or 50, got: {binder_contig}")
        self.params["rfd_contigs"] = binder_contig + "/0 " + target_contig

    def get_hotspots(self) -> str | None:
        return self.params.get("hotspot_residues")

    def set_hotspots(self, hotspots: str | None):
        self.params["hotspot_residues"] = hotspots

    def get_refolding_design_type(self) -> str:
        return "binder"

    def get_refolding_native_pdb_path(self, contig_index: int) -> Optional[str]:
        # In binder design, we don't compare with fixed input motif, so return None
        return None
