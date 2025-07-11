import sys
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import rdmolops

def solve():
    """
    This function identifies molecules from SMILES, checks their Crippen logP,
    and for those with logP > 1, calculates and sums their Eccentric Connectivity Indices.
    The detailed calculation for each qualifying molecule is printed.
    """
    # The three molecules depicted in the image are identified by their SMILES strings.
    # Reactant 1: 3-((2,4,5-trimethoxybenzyl)oxy)-5-fluoropyridin-2-amine
    # Reactant 2: 6-chloro-6,7-dihydro-5H-pyrrolo[3,4-d]pyrimidine
    # Product: 1-((5-fluoropyridin-2-yl)amino)-6,7-dihydro-5H-pyrrolo[3,4-d]pyrimidine
    molecules = {
        "Reactant 1": "COc1cc(OC)c(OC)cc1OCc2c(N)ncc(F)c2",
        "Reactant 2": "ClN1CC2=C(C1)N=CN=C2",
        "Product": "Fc1cncc(NN2CC3=C(C2)N=CN=C3)c1"
    }

    total_eci_sum = 0
    
    print("Analyzing depicted molecules...\n")

    for name, smiles in molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error parsing SMILES for {name}: {smiles}\n")
            continue

        # Step 1: Calculate Crippen logP
        logp = Crippen.MolLogP(mol)
        
        print(f"--- Molecule: {name} ---")
        print(f"SMILES: {smiles}")
        print(f"Crippen logP = {logp:.4f}")

        # Step 2: Filter by logP > 1
        if logp > 1:
            print("Crippen logP > 1, calculating Eccentric Connectivity Index (ECI).")
            
            # Step 3a: Add explicit hydrogens for the calculation
            mol_h = Chem.AddHs(mol)
            num_atoms = mol_h.GetNumAtoms()
            
            # Step 3b: Calculate all-pairs shortest path distance matrix
            dist_matrix = rdmolops.GetDistanceMatrix(mol_h)
            
            # Step 3c: Get valences (degrees) for all atoms
            valences = [atom.GetDegree() for atom in mol_h.GetAtoms()]
            
            # Step 3d: Get eccentricities for all atoms
            eccentricities = [int(max(row)) for row in dist_matrix]

            # Step 3e: Calculate ECI and build the equation string for output
            eci = 0
            equation_parts = []
            for i in range(num_atoms):
                v = valences[i]
                e = eccentricities[i]
                term = v * e
                eci += term
                equation_parts.append(f"{v} * {e}")

            # Print the detailed calculation equation as requested
            equation_str = " + ".join(equation_parts)
            print(f"ECI Calculation (Sum of Valence * Eccentricity for {num_atoms} atoms):")
            print(f"ECI = {equation_str}")
            print(f"ECI = {int(eci)}")
            
            total_eci_sum += eci
        else:
            print("Crippen logP <= 1, excluding from sum.")

        print("-" * (len(name) + 14) + "\n")

    # Step 4: Print the final sum
    print("====================================================================")
    print(f"Final Sum of Eccentric Connectivity Indices for molecules with logP > 1: {int(total_eci_sum)}")
    print("====================================================================")

solve()
<<<1814>>>