import sys
import subprocess

# Function to install a package using pip
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Try to import rdkit, and install it if it's not found
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen
    import numpy as np
except ImportError:
    print("RDKit not found. Installing it now...")
    install("rdkit-pypi")
    print("RDKit has been installed. Please run the script again.")
    sys.exit()

from rdkit import Chem
from rdkit.Chem import Crippen
import numpy as np

def calculate_eci_and_print(mols_info):
    """
    Calculates logP, filters molecules, computes ECI, and prints the full equation.
    """
    # Step 2: Calculate Crippen logP and filter molecules
    print("--- Step 1: Identifying molecules with Crippen logP > 1 ---")
    filtered_mols_info = []
    for name, smiles in mols_info:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: Could not parse SMILES for {name}: {smiles}")
            continue
        
        logp = Crippen.MolLogP(mol)
        print(f"Molecule '{name}': SMILES = {smiles}, logP = {logp:.3f}")
        
        if logp > 1:
            filtered_mols_info.append((name, mol))
            print(f"-> {name} has logP > 1 and is included in the calculation.")
        else:
            print(f"-> {name} has logP <= 1 and is excluded.")

    if not filtered_mols_info:
        print("\nNo molecules with Crippen logP > 1 found.")
        return

    # Step 3 & 4: Calculate Sum of Eccentric Connectivity Indices
    print("\n--- Step 2: Calculating Sum of Eccentric Connectivity Indices ---")
    total_eci_sum = 0
    all_terms_strings = []

    for name, mol in filtered_mols_info:
        # Add explicit hydrogens, as required by the problem
        mol_h = Chem.AddHs(mol)
        
        # Get the all-pairs shortest path distance matrix
        dist_matrix = Chem.GetDistanceMatrix(mol_h)
        
        num_atoms = mol_h.GetNumAtoms()
        mol_eci = 0
        
        for i in range(num_atoms):
            atom = mol_h.GetAtomWithIdx(i)
            # v_i: degree of the atom (including bonds to hydrogens)
            degree = atom.GetDegree()
            
            # e_i: eccentricity of the atom (max shortest path to any other atom)
            eccentricity = int(np.max(dist_matrix[i]))
            
            term = degree * eccentricity
            mol_eci += term
            all_terms_strings.append(f"{degree}*{eccentricity}")
        
        total_eci_sum += mol_eci

    # Step 5: Print the results as requested
    print("\nThe final calculation is the sum of (degree * eccentricity) for each atom")
    print("in all molecules that satisfy the logP > 1 condition.")
    print("\nFinal Equation:")
    
    # Print the equation in a readable format (10 terms per line)
    for i in range(0, len(all_terms_strings), 10):
        print(" + ".join(all_terms_strings[i:i+10]))
        
    print(f"\nSum of Eccentric Connectivity Indices = {total_eci_sum}")
    return total_eci_sum


if __name__ == "__main__":
    # Step 1: Define SMILES strings for the three molecules depicted in the image.
    # Reactant 1, Reactant 2, and Product
    molecules_to_analyze = [
        ("Reactant 1", "COc1cc(COC2=C(N)N=CC(F)=C2)cc(OC)c1OC"),
        ("Reactant 2", "c1ncc2n(Cl)cCC2n1"),
        ("Product", "Nc1cc(F)c(cn1)N1CCc2nccn2C1")
    ]
    
    final_sum = calculate_eci_and_print(molecules_to_analyze)
    
    # The final answer is submitted in the required format
    if final_sum is not None:
        print(f"\n<<< {final_sum} >>>")
