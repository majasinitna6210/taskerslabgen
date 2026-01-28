from pathlib import Path

from ase.build import surface
from ase.io import write


def build_cut_slabs(
    bulk_atoms,
    miller,
    layer_thickness_list,
    zbot,
    ztop,
    L,
    vacuum=10.0,
    out_dir=".",
    filename_template="surf_bulk_layers_{layer_thickness}.xyz",
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    for layer_thickness in layer_thickness_list:
        surf_bulk_n = surface(bulk_atoms, miller, layers=layer_thickness + 2, vacuum=0.0)
        zmin = zbot
        zmax = ztop + layer_thickness * L
        mask = [(zmin <= atom.position[2] <= zmax) for atom in surf_bulk_n]
        slab = surf_bulk_n[mask]
        slab.center(vacuum=vacuum, axis=2)

        out_path = out_dir / filename_template.format(layer_thickness=layer_thickness)
        write(out_path.as_posix(), slab)
        outputs.append(out_path.as_posix())

    return outputs
