# The user might need to install the RDKit library.
# This can be done by running the following command in the terminal:
# pip install rdkit

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
except ImportError:
    print("Error: The RDKit library is not installed.")
    print("Please install it using: pip install rdkit")
    # Provide pre-calculated results as a fallback
    print("\n--- Pre-calculated Information for Compound A ---")
    COMPOUND_A_NAME = "2-(3-hydroxypyridin-2-yl)-2-(phenylamino)acetonitrile"
    COMPOUND_A_FORMULA = "C13H11N3O"
    COMPOUND_A_MW = 225.25
    print(f"IUPAC Name: {COMPOUND_A_NAME}")
    print(f"Molecular Formula: {COMPOUND_A_FORMULA}")
    print(f"Number of Carbon atoms: 13")
    print(f"Number of Hydrogen atoms: 11")
    print(f"Number of Nitrogen atoms: 3")
    print(f"Number of Oxygen atoms: 1")
    print(f"Molecular Weight: {COMPOUND_A_MW:.3f} g/mol")

else:
    # The structure of Compound A, 2-(3-hydroxypyridin-2-yl)-2-(phenylamino)acetonitrile,
    # can be represented by the following SMILES string.
    smiles_A = "N#CC(Nc1ccccc1)c2c(O)cccn2"

    # Create a molecule object from the SMILES string
    mol_A = Chem.MolFromSmiles(smiles_A)

    if mol_A:
        # --- Calculate Properties ---
        # Molecular Formula
        formula = rdMolDescriptors.CalcMolFormula(mol_A)
        # Molecular Weight
        mol_weight = Descriptors.MolWt(mol_A)
        # IUPAC Name (hardcoded as RDKit does not have a free IUPAC name generator)
        iupac_name = "2-(3-hydroxypyridin-2-yl)-2-(phenylamino)acetonitrile"
        
        # --- Print Results ---
        print("--- Analysis of Final Product (Compound A) ---")
        print(f"IUPAC Name: {iupac_name}")
        print(f"SMILES String: {smiles_A}")
        print(f"Molecular Formula: {formula}")
        
        # Extract and print the number of each atom from the formula
        atom_counts = Chem.rdMolDescriptors.GetMolFormula(mol_A, separateInorg=False)
        # This gives a string like 'C13H11N3O', need to parse it if we want to print each number.
        # A simpler way is to iterate through atoms.
        element_counts = {}
        for atom in mol_A.GetAtoms():
            symbol = atom.GetSymbol()
            element_counts[symbol] = element_counts.get(symbol, 0) + 1
        
        print("\nAtom Counts:")
        print(f"Number of Carbon atoms: {element_counts.get('C', 0)}")
        print(f"Number of Hydrogen atoms: {element_counts.get('H', 0)}")
        print(f"Number of Nitrogen atoms: {element_counts.get('N', 0)}")
        print(f"Number of Oxygen atoms: {element_counts.get('O', 0)}")

        print(f"\nMolecular Weight: {mol_weight:.3f} g/mol")

    else:
        print("Error: Could not parse the SMILES string for Compound A.")
