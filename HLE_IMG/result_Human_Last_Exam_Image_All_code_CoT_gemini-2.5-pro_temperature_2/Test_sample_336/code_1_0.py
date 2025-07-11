import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen
import numpy as np

def calculate_eccentric_connectivity_index(mol_with_hs):
    """
    Calculates the Eccentric Connectivity Index for a molecule, including hydrogens.
    ECI = sum(degree(v) * eccentricity(v)) for all atoms v.
    """
    # Get the graph distance matrix (shortest path lengths between all atom pairs)
    dist_matrix = Chem.GetDistanceMatrix(mol_with_hs)
    num_atoms = mol_with_hs.GetNumAtoms()
    
    eci = 0
    eci_terms_list = []
    
    # Iterate over each atom (vertex) i
    for i in range(num_atoms):
        atom = mol_with_hs.GetAtomWithIdx(i)
        
        # Get the degree of the atom
        degree = atom.GetDegree()
        
        # Get the eccentricity of the atom (max value in its row of the distance matrix)
        eccentricity = int(max(dist_matrix[i]))
        
        # Calculate the term for this atom and add to sum
        term = degree * eccentricity
        eci += term
        eci_terms_list.append(f"({degree} * {eccentricity})")
        
    return eci, eci_terms_list

def solve_task():
    """
    Main function to solve the user's task.
    1. Define molecules by SMILES.
    2. Filter molecules by Crippen logP > 1.
    3. Calculate ECI for filtered molecules.
    4. Sum ECIs and print the result.
    """
    # The SMILES strings for the depicted molecules.
    # Note: Structures are interpreted from the image. Reactant 1 is simplified
    # to a representative structure due to ambiguity.
    smiles_definitions = {
        "Reactant 1": "c1ccc(cc1)COC1=C(N)C=C(F)C=N1",
        "Reactant 2": "ClN1CC2=C(C1)N=CC=N2",
        "Product": "NC1=C(N2CC3=C(C2)N=CC=N3)C=C(F)C=N1"
    }

    total_eci = 0
    all_eci_terms = []

    print("Analyzing depicted molecules...\n")
    for name, smiles in smiles_definitions.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Error parsing SMILES for {name}")
            continue

        # Add hydrogens, as required for the calculation
        mol_with_hs = Chem.AddHs(mol)

        # Calculate Crippen logP
        logp = Crippen.MolLogP(mol_with_hs)
        print(f"Molecule: {name}")
        print(f"SMILES: {smiles}")
        print(f"Crippen logP: {logp:.4f}")

        # Check if logP > 1
        if logp > 1:
            print("  -> logP > 1. Including in ECI calculation.\n")
            eci, eci_terms = calculate_eccentric_connectivity_index(mol_with_hs)
            total_eci += eci
            all_eci_terms.extend(eci_terms)
        else:
            print("  -> logP <= 1. Skipping molecule.\n")

    print("---")
    print("Final Calculation: Sum of Eccentric Connectivity Indices")
    if not all_eci_terms:
        print("No molecules found with Crippen logP > 1.")
        final_sum = 0
    else:
        # We need to print each term as number, not (deg*ecc)
        equation_numbers = [str(eval(term)) for term in all_eci_terms]
        equation_str = " + ".join(equation_numbers)
        print(f"Sum = {equation_str}")
        final_sum = sum([int(n) for n in equation_numbers])

    print(f"\nTotal Sum = {final_sum}")
    # Final answer block
    print(f"\n<<<{final_sum}>>>")

# Run the solver
solve_task()