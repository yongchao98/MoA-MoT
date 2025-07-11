import sys
import subprocess

# Install rdkit if it's not already installed
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdmolops
except ImportError:
    print("RDKit not found. Installing...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdmolops

def calculate_eci(mol):
    """
    Calculates the Eccentric Connectivity Index for a molecule.
    The molecule object should have explicit hydrogens added.
    """
    # Get the distance matrix (shortest path between all pairs of atoms)
    dist_matrix = rdmolops.GetDistanceMatrix(mol)
    
    # Get atom degrees
    degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
    
    # Get eccentricities (max value in each row of the distance matrix)
    eccentricities = [int(max(row)) for row in dist_matrix]
    
    # Calculate ECI and its individual terms
    eci_terms_str = []
    eci_value = 0
    for i in range(mol.GetNumAtoms()):
        term_value = degrees[i] * eccentricities[i]
        eci_terms_str.append(f"{degrees[i]}*{eccentricities[i]}")
        eci_value += term_value
        
    return eci_value, eci_terms_str

def main():
    """
    Main function to perform the requested calculations.
    """
    # SMILES strings for the three depicted molecules
    smiles_list = [
        "COc1cc(OC)c(c(OC)c1)COc1cncc(N)c1F",  # Molecule A (Reactant 1)
        "c1ncc2c(n1)CCN(Cl)C2",                # Molecule B (Reactant 2)
        "Fc1cc(NC2=NC=NC3=C2CCNC3)ccn1"        # Molecule C (Product)
    ]
    
    molecule_names = ["Molecule A (Reactant 1)", "Molecule B (Reactant 2)", "Molecule C (Product)"]
    
    molecules_to_process = []
    
    print("--- Analyzing Depicted Molecules ---")
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Could not parse SMILES for {molecule_names[i]}: {smiles}")
            continue
        
        # Add hydrogens for calculations
        mol_h = Chem.AddHs(mol)
        
        # Calculate Crippen logP
        logp = Descriptors.MolLogP(mol_h)
        
        print(f"\n{molecule_names[i]}:")
        print(f"  SMILES: {smiles}")
        print(f"  Crippen logP: {logp:.4f}")
        
        if logp > 1:
            print("  Result: logP > 1, molecule is INCLUDED in the sum.")
            molecules_to_process.append((molecule_names[i], mol_h))
        else:
            print("  Result: logP <= 1, molecule is EXCLUDED from the sum.")

    total_eci = 0
    
    print("\n--- Calculating Eccentric Connectivity Indices ---")
    
    if not molecules_to_process:
        print("No molecules with logP > 1 found.")
    else:
        for name, mol in molecules_to_process:
            eci, eci_terms = calculate_eci(mol)
            total_eci += eci
            
            print(f"\n{name}:")
            print(f"  Number of atoms (including H): {mol.GetNumAtoms()}")
            print(f"  Eccentric Connectivity Index (ECI) Calculation:")
            # To avoid extremely long output, we show a summary
            if len(eci_terms) > 30:
                 print(f"  {' + '.join(eci_terms[:15])} + ... + {' + '.join(eci_terms[-15:])} = {eci}")
            else:
                 print(f"  {' + '.join(eci_terms)} = {eci}")
            print(f"  ECI Value: {eci}")

    print("\n--- Final Result ---")
    print(f"The Sum of Eccentric Connectivity Indices for all molecules with logP > 1 is: {total_eci}")
    print(f"<<<{total_eci}>>>")

if __name__ == "__main__":
    main()