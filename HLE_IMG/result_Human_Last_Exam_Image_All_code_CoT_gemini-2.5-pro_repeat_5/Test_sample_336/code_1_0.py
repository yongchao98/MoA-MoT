import sys
import subprocess

# Install rdkit if it's not already installed
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen, AllChem
    import numpy as np
except ImportError:
    print("RDKit not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
    from rdkit import Chem
    from rdkit.Chem import Crippen, AllChem
    import numpy as np

def calculate_eci(mol, mol_name):
    """
    Calculates the Eccentric Connectivity Index (ECI) for a given RDKit molecule.
    ECI = sum(degree(i) * eccentricity(i)) for all atoms i.
    """
    # Add explicit hydrogens
    mol_h = Chem.AddHs(mol)
    
    # Get the graph distance matrix
    dist_matrix = Chem.GetDistanceMatrix(mol_h)
    
    eci = 0
    # Iterate over each atom to calculate its contribution to the ECI
    for i in range(mol_h.GetNumAtoms()):
        atom = mol_h.GetAtomWithIdx(i)
        degree = atom.GetDegree()
        # Eccentricity is the maximum value in the i-th row/column of the distance matrix
        eccentricity = int(np.max(dist_matrix[i]))
        eci += degree * eccentricity
        
    print(f"ECI for {mol_name}: {eci}")
    return eci

def solve():
    """
    Main function to solve the problem.
    """
    # SMILES strings for the three molecules depicted in the image
    smiles_list = [
        'Nc1nc(F)cc(OCc2cc(OC)c(OC)cc2OC)c1',  # Molecule 1 (left reactant)
        'c1ncc2c(n1)C[N](Cl)C2',                 # Molecule 2 (right reactant)
        'Nc1nc(F)c(N2CC3=C(C2)N=CN=C3)cc1'     # Molecule 3 (product)
    ]
    
    molecule_names = [
        "Molecule 1",
        "Molecule 2",
        "Molecule 3"
    ]

    print("Step 1: Calculate Crippen logP for each molecule.")
    
    molecules_to_process = []
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Error parsing SMILES for {molecule_names[i]}: {smiles}")
            continue
        
        logp = Crippen.MolLogP(mol)
        print(f"{molecule_names[i]} (SMILES: {smiles}) has a logP of {logp:.3f}.")
        
        if logp > 1:
            print(f"  -> logP > 1. Including this molecule for ECI calculation.\n")
            molecules_to_process.append((molecule_names[i], mol))
        else:
            print(f"  -> logP <= 1. Excluding this molecule.\n")
    
    print("Step 2: Calculate Eccentric Connectivity Index (ECI) for the filtered molecules.")
    if not molecules_to_process:
        print("No molecules with logP > 1 were found.")
        total_eci = 0
        eci_values = []
    else:
        eci_values = [calculate_eci(mol, name) for name, mol in molecules_to_process]
        total_eci = sum(eci_values)

    print("\nStep 3: Calculate the Sum of Eccentric Connectivity Indices.")
    if len(eci_values) > 0:
        calculation_str = " + ".join(map(str, eci_values))
        print(f"The final Sum of Eccentric Connectivity Indices is: {calculation_str} = {total_eci}")
    else:
        print("The final Sum of Eccentric Connectivity Indices is: 0")
        
    return total_eci

if __name__ == "__main__":
    final_answer = solve()
    print(f"\n<<< {final_answer} >>>")
