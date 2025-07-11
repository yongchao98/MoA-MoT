import sys
import subprocess

# Install rdkit if it is not already installed
try:
    from rdkit import Chem
    from rdkit.Chem import Crippen
except ImportError:
    print("RDKit library not found. Installing it now...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
        from rdkit import Chem
        from rdkit.Chem import Crippen
        print("RDKit installed successfully.")
    except Exception as e:
        print(f"Failed to install RDKit. Please install it manually using 'pip install rdkit-pypi'. Error: {e}")
        sys.exit(1)

def calculate_eccentric_connectivity_index(mol):
    """
    Calculates the Eccentric Connectivity Index for a given RDKit molecule.
    The calculation includes hydrogen atoms.
    """
    # Add explicit hydrogens to the molecular graph for the calculation
    mol_with_hs = Chem.AddHs(mol)
    
    num_atoms = mol_with_hs.GetNumAtoms()
    if num_atoms == 0:
        return 0

    # Calculate the all-pairs shortest path distance matrix
    dist_matrix = Chem.GetDistanceMatrix(mol_with_hs)
    
    eci = 0
    # Iterate over each atom in the molecule (including hydrogens)
    for i in range(num_atoms):
        # Get the degree of the atom (number of connections)
        degree = mol_with_hs.GetAtomWithIdx(i).GetDegree()
        
        # Get the eccentricity of the atom (max distance to any other atom)
        if len(dist_matrix[i]) > 0:
            eccentricity = max(dist_matrix[i])
        else:
            eccentricity = 0
            
        # Add the product to the total ECI sum
        eci += degree * eccentricity
        
    return int(eci)

def solve():
    """
    Main function to identify molecules, filter by logP,
    and calculate the sum of Eccentric Connectivity Indices.
    """
    # SMILES strings for the molecules depicted in the image
    molecules = {
        "Reactant 1": "Nc1nccc(F)c1OCc1cc(OC)c(OC)c(OC)c1",
        "Reactant 2": "ClN1Cc2cncnc2C1",
        "Product": "Nc1nccc(F)c1N1Cc2cncnc2C1"
    }

    total_eci = 0
    eci_values = []

    print("Analyzing depicted molecules...\n")
    
    for name, smiles in molecules.items():
        print(f"--- Analyzing {name} ---")
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            print(f"Error: Could not parse SMILES string '{smiles}'")
            continue
        
        # Calculate Crippen logP
        logP = Crippen.MolLogP(mol)
        print(f"SMILES: {smiles}")
        print(f"Crippen logP: {logP:.2f}")

        # Check if the molecule meets the condition logP > 1
        if logP > 1:
            print("Condition (logP > 1) is MET.")
            # Calculate Eccentric Connectivity Index
            eci = calculate_eccentric_connectivity_index(mol)
            print(f"Eccentric Connectivity Index (ECI): {eci}")
            total_eci += eci
            eci_values.append(str(eci))
        else:
            print("Condition (logP > 1) is NOT MET. Skipping ECI calculation for this molecule.")
        
        print("-" * (len(name) + 16) + "\n")

    print("=========================================================")
    print("Final Result:")
    print("The sum of Eccentric Connectivity Indices for molecules with Crippen logP > 1 is:")
    
    if not eci_values:
        print("No molecules met the required condition.")
        print("Final Sum = 0")
        total_eci = 0
    else:
        # Format the final equation as requested
        equation_str = " + ".join(eci_values)
        print(f"{equation_str} = {total_eci}")
    
    return total_eci

if __name__ == "__main__":
    final_answer = solve()
    # The final answer is wrapped for the system to parse.
    # print(f"\n<<< {final_answer} >>>") # This line will be executed if run as script
    # The following line is for direct output in this environment.
    print(f"\n<<<{final_answer}>>>")

# Calling the function to execute the plan.
# Note: RDKit installation may take a moment on first run.
solve()