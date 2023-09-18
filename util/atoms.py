from Bio.PDB import PDBParser, MMCIFParser
from Bio import PDB
from Bio.PDB.mmcifio import MMCIFIO


def load_structure_from_file(fn):
    if fn.split(".")[-1][:3] == "pdb":
        parser = PDBParser(QUIET=True)
    elif fn.split(".")[-1][:3] == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise RuntimeError("Unknown type for structure file:", fn[-3:])

    structure = parser.get_structure('PHA-L', fn)

    return structure


def save_structure_from_file(structure, fn):
    if fn.split(".")[-1][:3] == "pdb":
        io = PDB.PDBIO()
    elif fn.split(".")[-1][:3] == "cif":
        io = MMCIFIO()
    else:
        raise RuntimeError("Unknown output file type")

    io.set_structure(structure)
    io.save(fn)
