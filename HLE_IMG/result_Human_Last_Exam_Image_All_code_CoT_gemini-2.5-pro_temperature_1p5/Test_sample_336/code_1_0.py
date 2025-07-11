import sys
from io import StringIO

def install_and_import():
    """Installs and imports the RDKit library."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen
        return Chem, Crippen
    except ImportError:
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi"])
        from rdkit import Chem
        from rdkit.Chem import Crippen
        return Chem, Crippen

def calculate_eci(mol_with_hs):
    """
    Calculates the Eccentric Connectivity Index for a molecule.
    The molecule object must have explicit hydrogens added.
    ECI = sum(degree(i) * eccentricity(i)) for all atoms i.
    """
    dist_matrix = Chem.GetDistanceMatrix(mol_with_hs)
    num_atoms = mol_with_hs.GetNumAtoms()
    eci = 0
    for i in range(num_atoms):
        atom = mol_with_hs.GetAtomWithIdx(i)
        degree = atom.GetDegree()
        # Eccentricity is the maximum value in the atom's row/column of the distance matrix
        eccentricity = max(dist_matrix[i])
        eci += degree * eccentricity
    return int(eci)

def main():
    """
    Main function to identify molecules with logP > 1 and calculate the sum of their ECIs.
    """
    Chem, Crippen = install_and_import()

    # SMILES strings for the three molecules depicted in the image
    molecules = {
        "Reactant 1": "COc1cc(OC)c(OC)c(COc2c(N)nc(F)cc2)c1",
        "Reactant 2": "c1ncc2c(n1)CN(Cl)C2",
        "Product": "Nc1cc(N2Cc3c(CN2)ncc-3)c(F)cn1",
    }

    eci_values = []
    
    for name, smi in molecules.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            # This should not happen with valid SMILES
            continue

        # Calculate Crippen logP
        logp = Crippen.MolLogP(mol)
        
        # Check if the molecule meets the logP > 1 criterion
        if logp > 1:
            # Add explicit hydrogens for the ECI calculation
            mol_h = Chem.AddHs(mol)
            
            # Calculate the Eccentric Connectivity Index
            eci = calculate_eci(mol_h)
            eci_values.append(eci)

    # Calculate the sum and print the final equation
    total_eci = sum(eci_values)
    equation_str = " + ".join(map(str, eci_values))
    print(f"The ECI values for molecules with logP > 1 are: {', '.join(map(str, eci_values))}")
    print("The final sum is calculated as follows:")
    print(f"{equation_str} = {total_eci}")
    
    # Store final answer for parsing
    # Use StringIO to capture the output and format it for the final response
    output_catcher = StringIO()
    sys.stdout = output_catcher
    print(total_eci)
    sys.stdout = sys.__stdout__ # Restore stdout
    final_answer = output_catcher.getvalue().strip()
    print(f"\n<<<${final_answer}$>>>")


if __name__ == "__main__":
    main()