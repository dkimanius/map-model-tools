from Bio.PDB import PDBParser, MMCIFParser


def load_structure_from_file(fn):
    if fn.split(".")[-1][:3] == "pdb":
        parser = PDBParser(QUIET=True)
    elif fn.split(".")[-1][:3] == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise RuntimeError("Unknown type for structure file:", fn[-3:])

    structure = parser.get_structure('PHA-L', fn)

    return structure
