# To run this code, you need to have the RDKit library installed.
# You can install it using pip: pip install rdkit
# Or using conda: conda install -c conda-forge rdkit

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    import re
except ImportError:
    print("RDKit library not found.")
    print("Please install it to run this code:")
    print("pip install rdkit")
    exit()

def identify_compound_c():
    """
    This function identifies and describes the final product, Compound C.
    """
    # The structure of Compound C is deduced as tris(4-diethylamino-2,6-dihydroxyphenyl)methylium.
    # We represent it using a SMILES string.
    # [C+] is the central cationic carbon.
    # It's attached to three identical (c1...) rings.
    # c1c(O)cc(N(CC)CC)cc1O describes a 4-diethylamino-2,6-dihydroxyphenyl group.
    
    smiles_c = "CCN(CC)c1cc(O)c(c(O)c1)[C+](c1c(O)cc(N(CC)CC)cc1O)c1c(O)cc(N(CC)CC)cc1O"
    
    # Create a molecule object from the SMILES string
    mol_c = Chem.MolFromSmiles(smiles_c)
    
    if mol_c:
        # Add explicit hydrogens for accurate formula calculation
        mol_c_with_h = Chem.AddHs(mol_c)
        
        # Calculate the molecular formula
        formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol_c_with_h)
        
        # Remove the charge from the formula string for parsing
        formula_neutral = formula.replace("+", "").replace("-", "")

        print("--- Analysis of Final Product (Compound C) ---")
        print("\nChemical Name: tris(4-diethylamino-2,6-dihydroxyphenyl)methylium")
        print(f"SMILES String: {smiles_c}")
        print(f"\nMolecular Formula (for the cation): {formula}")
        
        print("\nAtom count in the final molecule (cation):")
        
        # Use regex to find all element-number pairs
        atom_counts = re.findall('([A-Z][a-z]*)(\d*)', formula_neutral)
        
        for element, count in atom_counts:
            # If count is an empty string, it means 1
            num = int(count) if count else 1
            print(f"- {element} (Carbon): {num}" if element == 'C' else f"- {element}: {num}")

    else:
        print("Error: Could not parse the SMILES string for Compound C.")

if __name__ == "__main__":
    identify_compound_c()
