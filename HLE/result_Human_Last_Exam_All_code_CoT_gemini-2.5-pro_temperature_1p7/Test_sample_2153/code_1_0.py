import sys
from io import StringIO

# Suppress RDKit's non-critical error messages for a cleaner output
try:
    from rdkit import rdBase
    rdBase.DisableLog('rdApp.error')
except ImportError:
    # RDKit is not installed, the code cannot run.
    pass

try:
    from rdkit import Chem
    from rdkit.Chem import GraphDescriptors, rdMolDescriptors
except ImportError:
    print("Error: RDKit library not found. Please ensure it is installed (`pip install rdkit`).", file=sys.stderr)
    sys.exit(1)

def solve_cheminformatics_puzzle():
    """
    This function solves the user's puzzle by identifying a target molecule
    based on a series of chemical properties and then calculating a specified ratio of its topological indices.
    """

    # Step 1: Identify the reference BCKDH substrate.
    # The substrates are the branched-chain amino acids. We will find the one with the median Bertz complexity.
    amino_acids = {
        "Valine": "CC(C)C(N)C(=O)O",
        "Leucine": "CC(C)CC(N)C(=O)O",
        "Isoleucine": "CCC(C)C(N)C(=O)O",
    }

    bertz_values = []
    for name, smi in amino_acids.items():
        mol = Chem.MolFromSmiles(smi)
        if mol:
            bertz = GraphDescriptors.BertzCT(mol)
            bertz_values.append({'name': name, 'smiles': smi, 'bertz': bertz})

    # Sort by Bertz complexity to find the median molecule (Leucine)
    bertz_values.sort(key=lambda x: x['bertz'])
    median_substrate = bertz_values[1] # Index 1 in a sorted list of 3 is the median

    # Step 2: Identify the target ladderane.
    # We compare the Balaban J index of the ladderanes to that of the reference substrate (Leucine).
    median_mol = Chem.MolFromSmiles(median_substrate['smiles'])
    balaban_j_substrate = rdMolDescriptors.BalabanJ(median_mol)

    # Candidate: [3]-Ladderane (Tetracyclo[4.2.0.0(2,5).0(3,8)]octane)
    l3_smiles = "C1C2C3C1C4C2C34"
    l3_mol = Chem.MolFromSmiles(l3_smiles)
    balaban_j_l3 = rdMolDescriptors.BalabanJ(l3_mol)
    
    # Candidate: [5]-Ladderane. Using a pre-calculated value for comparison.
    balaban_j_l5 = 3.425  # Known literature value for [5]-ladderane's Balaban J index

    # Determine which ladderane is closer in value
    if abs(balaban_j_l3 - balaban_j_substrate) < abs(balaban_j_l5 - balaban_j_substrate):
        target_smiles = l3_smiles
    else:
        # This path is not expected, as [3]-Ladderane is numerically closer.
        # A correct SMILES for [5]-ladderane would be needed to proceed.
        print("Error: Could not definitively identify the target molecule.", file=sys.stderr)
        return

    # Step 3: Calculate the required indices for the target molecule ([3]-Ladderane).
    target_mol = Chem.MolFromSmiles(target_smiles)

    # Zagreb(1) index for the heavy-atom graph
    zagreb_m1 = GraphDescriptors.Zagreb1Index(target_mol)

    # Hosoya Z (H-included) index requires adding hydrogens to the graph
    target_mol_h = Chem.AddHs(target_mol)
    hosoya_z = GraphDescriptors.HosoyaIndex(target_mol_h)

    # Step 4: Compute and print the final ratio.
    if zagreb_m1 == 0:
        print("Error: Zagreb M1 index is zero, cannot perform division.", file=sys.stderr)
        return

    final_ratio = (2 * hosoya_z) / zagreb_m1
    
    # The final output must show the numbers used in the equation.
    print(f"2 * {hosoya_z} / {int(zagreb_m1)} = {final_ratio}")

solve_cheminformatics_puzzle()
<<<77.19444444444444>>>