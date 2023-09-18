import argparse

import numpy as np
import torch
from Bio import PDB

from util.atoms import load_structure_from_file
from util.grid import load_mrc, sample_from_map

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('mrc', type=str, help="MRC file")
    parser.add_argument('input', type=str, help="Input structure file")
    parser.add_argument('output', type=str, help="output structure file")
    args = parser.parse_args()

    grid_values, voxel_size, global_origin = load_mrc(args.mrc)
    grid_values = torch.from_numpy(grid_values).float()

    structure = load_structure_from_file(args.input)
    if len(structure) > 1:
        print(f"WARNING: {len(structure)} structures found in model file: {args.input}")
    coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_vector().get_array())
    coordinates = torch.from_numpy(np.array(coordinates)).float()

    sampled_values = sample_from_map(grid_values, global_origin, voxel_size, coordinates)

    print(f"Min={sampled_values.min()}")
    print(f"max={sampled_values.max()}")
    print(f"mean={sampled_values.mean()}")

    sampled_values = sampled_values.cpu().detach().numpy()

    # Iterate through atoms and set B-factors
    i = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    # You can set the B-factor value as desired
                    atom.set_bfactor(sampled_values[i])
                    i += 1

    # Save the modified structure to a new PDB file
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(args.output)
