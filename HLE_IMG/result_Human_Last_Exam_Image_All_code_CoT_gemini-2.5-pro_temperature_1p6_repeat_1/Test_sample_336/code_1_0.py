import rdkit
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import GraphDescriptors

def solve_chemistry_problem():
    """
    This script identifies molecules with Crippen logP > 1 from a given set,
    calculates their Eccentric Connectivity Index (ECI), and prints the sum.
    The ECI calculation includes hydrogen atoms.
    """
    # SMILES strings for the three molecules based on visual interpretation of the reaction
    # Reactant 1: 3-((2,3,4-trimethoxybenzyl)oxy)-5-fluoropyridin-2-amine
    smiles_r1 = "COC1=C(OC)C(OC)=CC=C1COC2=C(N)C(=CN=C2)F"
    # Reactant 2: 6-chloro-6,7-dihydro-5H-pyrrolo[3,4-b]pyrazine
    smiles_r2 = "c1cnc2c(n1)CN(Cl)C2"
    # Product: N-(6,7-dihydro-5H-pyrrolo[3,4-b]pyrazin-2-yl)-5-fluoropyridin-2-amine
    smiles_p = "Fc1cc(Nc2cnc3c(n2)NCC3)cn1"

    molecules = {
        "Reactant 1": smiles_r1,
        "Reactant 2": smiles_r2,
        "Product": smiles_p
    }

    qualifying_molecules = []

    print("--- Analysis of Depicted Molecules ---")
    for name, smiles in molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Could not parse SMILES for {name}")
            continue
        
        mol_h = Chem.AddHs(mol)
        logp = Crippen.MolLogP(mol_h)
        print(f"{name} (SMILES: {smiles}) has Crippen logP = {logp:.3f}")
        
        if logp > 1:
            print(f"-> This molecule qualifies (logP > 1).")
            qualifying_molecules.append(mol_h)
        else:
            print(f"-> This molecule does not qualify (logP <= 1).")
    
    print("\n--- Calculating Sum of Eccentric Connectivity Indices ---")

    total_eci = 0
    all_terms_for_print = []

    if not qualifying_molecules:
        print("No molecules with logP > 1 found.")
    else:
        for mol_h in qualifying_molecules:
            dist_matrix = Chem.GetDistanceMatrix(mol_h)
            num_atoms = mol_h.GetNumAtoms()
            for i in range(num_atoms):
                atom = mol_h.GetAtomWithIdx(i)
                degree = atom.GetDegree()
                # Eccentricity for atom i is the max value in its row/column in the distance matrix
                eccentricity = int(max(dist_matrix[i]))
                term_value = degree * eccentricity
                total_eci += term_value
                all_terms_for_print.append(f"{degree}*{eccentricity}")

        equation_str = " + ".join(all_terms_for_print)
        print("The final sum is calculated from the following terms (degree * eccentricity for each atom):")
        print(f"{equation_str} = {total_eci}")

    # Final answer in the required format
    print(f"\nFinal calculated sum: {total_eci}")
    print(f"<<<{total_eci}>>>")

solve_chemistry_problem()