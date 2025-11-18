plugin = dict(
    pages = {
        "🪩 ProteinDJ": [
            dict(page="ovo_proteindj.pages.binder_denovo", title="🖇️ ProteinDJ De Novo Binder"),
            dict(page="ovo_proteindj.pages.monomer_motifscaff", title="🏗 ProteinDJ Scaffold Design"),
        ],
    },
    descriptors="ovo_proteindj.descriptors_proteindj",
    modules = [
        "ovo_proteindj.models_proteindj"
    ]
)
