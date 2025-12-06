import os
import re
from dataclasses import dataclass, field
from typing import Callable, Literal

from ovo import config
from ovo.core.database import DesignWorkflow, WorkflowTypes, WorkflowParams, Base
from ovo.core.utils.residue_selection import from_segments_to_hotspots, from_hotspots_to_segments, from_contig_to_residues

plugin_config = config.plugins.get(__package__, {})
default_params = plugin_config.get("default_params", {})

@WorkflowTypes.register("ProteinDJ")
@dataclass
class ProteinDJDesignWorkflow(DesignWorkflow):

    @dataclass
    class Params(WorkflowParams):
        input_pdb: str = None
        rfd_contigs: str = "[]"
        design_mode: str = None
        seq_method: str = "mpnn"
        pred_method: str = "af2"
        out_dir: str = "output"
        num_designs: int = 10
        seqs_per_design: int = 8
        rfd_extra_config: str = ""
        # Number of available GPU machines to be used in parallel - determines batch size
        gpus: int = 1
        # Parameters recommended to be set through ovo config
        rfd_models: str = default_params.get("rfd_models").format(config=config)
        af2_models: str = default_params.get("af2_models").format(config=config)
        boltz_models: str = default_params.get("boltz_models").format(config=config)
        cpus: int = default_params.get("cpus")
        cpus_per_gpu: int = default_params.get("cpus_per_gpu")
        memory_gpu: str = default_params.get("memory_gpu")
        memory_cpu: str = default_params.get("memory_cpu")
        gpu_queue: str = default_params.get("gpu_queue")
        gpu_model: str = default_params.get("gpu_model")

        # FIXME !!!
        # FIXME remove this - required in schema but should be null by default
        # FIXME !!!
        design_length: str = "123-456"
        rfd_partial_diffusion_timesteps: int = 20
        rfd_scaffold_dir: str = "UNUSED"
        # FIXME !!!

        def validate(self):
            if not self.input_pdb:
                raise ValueError("input_pdb is required")
            if not self.rfd_contigs or self.rfd_contigs == "[]":
                raise ValueError("rfd_contigs is required and cannot be empty")
            if self.num_designs is None or self.num_designs < 1:
                raise ValueError("num_designs must be a positive integer")
            if self.seqs_per_design is None or self.seqs_per_design < 1:
                raise ValueError("seqs_per_design must be a positive integer")

    params: Params = field(default_factory=Params, metadata=dict(tool_name="ProteinDJ"))
    preview_job_id: str = None

    def get_pipeline_name(self) -> str:
        return "https://github.com/PapenfussLab/proteindj@main"

    def prepare_params(self, workdir: str) -> dict:
        from ovo import storage
        params = self.params.to_dict()
        params["input_pdb"] = storage.prepare_workflow_input(params["input_pdb"], workdir)
        return params

    def process_results(self, job: "DesignJob", callback: Callable = None) -> list[Base]:
        """Process results of a successful workflow - download files from workdir, save DesignJob, Pool and Designs"""
        from ovo_proteindj.logic import process_workflow_results
        return process_workflow_results(job=job, callback=callback)

    def get_input_name(self) -> str | None:
        return os.path.basename(self.params.input_pdb.rsplit(".", 1)[0]) if self.params.input_pdb else None

    def get_input_pdb_path(self) -> str | None:
        return self.params.input_pdb

    def get_selected_segments(self) -> list[str] | None:
        raise NotImplementedError()

    def set_selected_segments(self, segments: list[str]):
        raise NotImplementedError()

    def get_contig(self, contig_index: int = 0) -> list[str] | Literal[""]:
        contigs = self.params.rfd_contigs.strip("[]").split(",") if self.params.rfd_contigs else []
        return contigs[contig_index] if contigs else ""

    def set_contig(self, contig: str):
        self.params.rfd_contigs = f"[{contig.strip('[]')}]" if contig else "[]"

    def get_cyclic_offset(self) -> bool:
        return False

    def get_hotspots(self) -> str | None:
        return None



@WorkflowTypes.register("ProteinDJ monomer_denovo")
@dataclass
class ProteinDJMonomerMotifScaffDesignWorkflow(ProteinDJDesignWorkflow):

    @dataclass
    class Params(ProteinDJDesignWorkflow.Params):
        design_mode: str = "monomer_motifscaff"
        seq_method: str = "fampnn"

    params: Params = field(default_factory=Params, metadata=dict(tool_name="ProteinDJ"))

    def get_inpaint_seq(self) -> str | None:
        # FIXME parse from params.rfd_extra_config
        raise NotImplementedError()

    def set_inpaint_seq(self):
        # FIXME add to params.rfd_extra_config
        raise NotImplementedError()

    def get_selected_segments(self, inpainting=False) -> list[str] | None:
        if inpainting:
            raise NotImplementedError("Inpainting not implemented yet")
        contig = self.get_contig()
        return [segment for subcontig in contig.split() for segment in subcontig.split("/") if segment and segment[0].isalpha()]

    def set_selected_segments(self, segments: list[str]):
        # TODO preserve generated segments, see ovo logic for re-generating contig
        # contig = self.get_contig()
        # generated_segments = [segment for subcontig in contig.split() for segment in subcontig.split("/") if segment and segment[0].isnumeric()]
        # if generated_segments:
        #     see ovo logic for re-generating contig
        # else:
        self.set_contig("/".join(segments))


@WorkflowTypes.register("ProteinDJ binder_denovo")
@dataclass
class ProteinDJBinderDeNovoDesignWorkflow(ProteinDJDesignWorkflow):

    @dataclass
    class Params(ProteinDJDesignWorkflow.Params):
        hotspot_residues: str = None
        design_mode: str = "binder_denovo"


    params: Params = field(default_factory=Params, metadata=dict(tool_name="ProteinDJ"))

    preview_job_id: str = None

    # TODO adjust this for scaffold vs binder
    def get_selected_segments(self) -> list[str] | None:
        return from_hotspots_to_segments(self.params.hotspot_residues)

    def set_selected_segments(self, segments: list[str]) -> None:
        # Convert from segments to hotspots
        residues = from_segments_to_hotspots(segments)
        self.set_hotspots(residues)

    def get_target_chain(self) -> str | None:
        if not self.params.rfd_contigs or self.params.rfd_contigs == "[]":
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
        self.params.rfd_contigs = binder_contig + " " + target_contig

    def set_binder_contig(self, binder_contig: str):
        target_contig = self.get_target_contig()
        assert target_contig, "Target contig must be set before setting binder contig"
        if not all(re.fullmatch(r"^\d+(-\d+)?$", segment) for segment in binder_contig.split("/")):
            raise ValueError(f"Expected binder contig to be in format 50-100 or 50, got: {binder_contig}")
        self.params.rfd_contigs = binder_contig + "/0 " + target_contig

    def get_hotspots(self) -> str | None:
        return self.params.hotspot_residues

    def set_hotspots(self, hotspots: str | None):
        self.params.hotspot_residues = hotspots
