import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def get_eci_terms(smiles_string):
    """
    Calculates the terms (degree * eccentricity) for the ECI of a molecule.
    Each term corresponds to one atom (including hydrogens).

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        tuple: A tuple containing:
            - int: The total Eccentric Connectivity Index value.
            - list: A list of strings, where each string is a "degree * eccentricity" term.
    """
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        print(f"Warning: Could not parse SMILES '{smiles_string}'")
        return 0, []

    # The calculation must include hydrogens
    mol_with_hs = Chem.AddHs(mol)
    
    # The distance matrix stores the shortest path between any two atoms
    dist_matrix = Chem.GetDistanceMatrix(mol_with_hs)
    num_atoms = mol_with_hs.GetNumAtoms()

    eci = 0
    terms = []

    for i in range(num_atoms):
        atom = mol_with_hs.GetAtomWithIdx(i)
        degree = atom.GetDegree()
        # Eccentricity is the maximum value in the atom's row of the distance matrix
        eccentricity = int(np.max(dist_matrix[i]))
        
        term_value = degree * eccentricity
        eci += term_value
        terms.append(f"{degree} * {eccentricity}")

    return eci, terms

def main():
    """
    Main function to execute the plan: identify, filter, calculate, and sum.
    """
    # SMILES strings for the three molecules identified from the image
    # Reactant 1: 3-((3,4,5-trimethoxybenzyl)oxy)-5-fluoropyridin-2-amine
    smiles_A = "Nc1nc(F)cc(OOCc2cc(OC)c(OC)c(OC)c2)c1"
    # Reactant 2: 7-chloro-6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidine
    smiles_B = "c1ncc2c(n1)N(Cl)CC2"
    # Product: 4-(6,7-dihydro-5H-pyrrolo[2,3-d]pyrimidin-7-yl)-5-fluoropyridin-2-amine
    smiles_C = "Nc1ncc(F)cc1N1CCc2ncncc12"

    all_molecules = [smiles_A, smiles_B, smiles_C]
    
    qualified_molecules_smiles = []
    print("Step 1: Filtering molecules with Crippen logP > 1")
    for i, smiles in enumerate(all_molecules):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            logP = Descriptors.MolLogP(mol)
            print(f"  - Molecule {i+1}: logP = {logP:.3f}", end="")
            if logP > 1:
                qualified_molecules_smiles.append(smiles)
                print(" -> Qualifies")
            else:
                print(" -> Does not qualify")
    
    print("\nStep 2: Calculating the Sum of Eccentric Connectivity Indices for qualified molecules.")

    grand_total_eci = 0
    all_terms = []

    for smiles in qualified_molecules_smiles:
        eci, terms = get_eci_terms(smiles)
        grand_total_eci += eci
        all_terms.extend(terms)
    
    print("\nThe final result is the sum of (degree * eccentricity) for every atom in every qualifying molecule.")
    print("\nFinal Equation:")
    # Printing the full equation as requested
    equation_string = " + ".join(all_terms)
    print(f"Sum = {equation_string}")
    
    print(f"\nTotal Sum of Eccentric Connectivity Indices = {grand_total_eci}")

if __name__ == "__main__":
    main()