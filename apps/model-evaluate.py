import argparse

import numpy as np
import torch
from Bio import PDB
from sklearn.neighbors import KDTree

from util.atoms import load_structure_from_file
from util.grid import load_mrc, sample_from_map

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('model', type=str, help="MRC file")
    parser.add_argument('ground_truth', type=str, help="Input structure file")
    parser.add_argument('--radius', type=float, default=3., help="Search max radius")
    args = parser.parse_args()

    structure = load_structure_from_file(args.model)
    if len(structure) > 1:
        print(f"WARNING: {len(structure)} structures found in model file: {args.model}")

    resnames = []
    coordinates = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resnames.append(residue.resname)
                coordinates.append(residue['CA'].get_vector().get_array())

    resnames = np.array(resnames)
    coordinates = np.array(coordinates)
    tree = KDTree(coordinates, leaf_size=2)

    coordinate_hit = 0.
    resname_hits = 0.
    total_count = 0.

    structure_target = load_structure_from_file(args.ground_truth)
    if len(structure_target) > 1:
        print(f"WARNING: {len(structure_target)} target structures found in model file: {args.ground_truth}")

    for model in structure_target:
        for chain in model:
            for residue in chain:
                if 'CA' not in residue:
                    continue
                coordinate = residue['CA'].get_vector().get_array()[None]

                resname = residue.resname
                total_count += 1

                hits = tree.query_radius(coordinate, r=args.radius, count_only=False)[0]

                if len(hits) > 0:
                    coordinate_hit += 1

                for i in hits:
                    if resnames[i] == resname:
                        resname_hits += 1

    print(f"total_count={int(total_count)}")
    print(f"coordinate_hit={int(coordinate_hit)}")
    print(f"resname_hits={int(resname_hits)}")
    print(f"sequence recovery={round(resname_hits/total_count, 3)}")
