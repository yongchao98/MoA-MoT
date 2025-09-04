# First, ensure you have rdkit installed:
# pip install rdkit-pypi
# pip install cirpy

import sys
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("RDKit not found. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)

try:
    import cirpy
except ImportError:
    print("cirpy not found. Please install it using 'pip install cirpy'")
    sys.exit(1)

def get_smiles_from_name(name):
    """Converts an IUPAC name to a SMILES string using the CIR service."""
    try:
        return cirpy.resolve(name, 'smiles')
    except Exception as e:
        print(f"Could not resolve name '{name}': {e}")
        return None

def check_answer():
    """
    Checks the correctness of the proposed answer for the chemistry problem.
    """
    # --- Step 1: Define molecules and verify their structures ---
    start_name = "3,3,6-trimethylhepta-1,5-dien-4-one"
    product_c_name = "5-hydroxy-3,3,6,6-tetramethylhept-1-en-4-one"

    start_smiles = get_smiles_from_name(start_name)
    product_c_smiles = get_smiles_from_name(product_c_name)

    if not all([start_smiles, product_c_smiles]):
        return "Failure: Could not convert one or more IUPAC names to structures. Cannot proceed."

    start_mol = Chem.MolFromSmiles(start_smiles)
    product_c_mol = Chem.MolFromSmiles(product_c_smiles)

    # --- Step 2: Analyze Reaction 1 (Epoxidation Site) ---
    # The starting material has two double bonds.
    # We expect epoxidation at the more electron-rich, more substituted C5=C6 bond.
    # C1=C2 is monosubstituted and conjugated (electron-poor).
    # C5=C6 is trisubstituted (electron-rich).
    
    # Find the double bonds by SMARTS pattern
    double_bonds = start_mol.GetSubstructMatches(Chem.MolFromSmarts('C=C'))
    if len(double_bonds) != 2:
        return f"Incorrect starting material structure: Expected 2 double bonds, found {len(double_bonds)}."

    # Identify the conjugated double bond C1=C2
    conjugated_db_pattern = Chem.MolFromSmarts('[CH2]=[CH]-C(=O)')
    if not start_mol.HasSubstructMatch(conjugated_db_pattern):
        return "Constraint Check Failed: The starting material structure does not match the expected α,β-unsaturated ketone system."

    # Identify the isolated, trisubstituted double bond C5=C6
    isolated_db_pattern = Chem.MolFromSmarts('C=C(C)C')
    if not start_mol.HasSubstructMatch(isolated_db_pattern):
        return "Constraint Check Failed: The starting material structure does not match the expected trisubstituted double bond."

    # Conclusion for Step 1: The presence of a more substituted, non-conjugated double bond
    # makes it the kinetically favored site for epoxidation. The logic is sound.
    
    # --- Step 3: Analyze Reaction 2 (Gilman Addition) ---
    
    # 3a: Check atom balance. The net reaction is epoxidation (+O) followed by
    # addition of a methyl group (CH3) from the Gilman reagent and a proton (H) from workup.
    # Total addition: O + CH4 = CH4O
    start_formula = CalcMolFormula(start_mol)
    product_c_formula = CalcMolFormula(product_c_mol)
    
    expected_product_formula = CalcMolFormula(Chem.MolFromSmiles(start_smiles + ".O.C")) # Add O and C
    # Manually adjust H count for CH4O addition
    start_h_count = Descriptors.HeavyAtomCount(start_mol) # This is not H count
    start_h_count = sum(1 for atom in start_mol.GetAtoms() if atom.GetAtomicNum() == 1)
    start_c_count = sum(1 for atom in start_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    start_o_count = sum(1 for atom in start_mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    prod_h_count = sum(1 for atom in product_c_mol.GetAtoms() if atom.GetAtomicNum() == 1)
    prod_c_count = sum(1 for atom in product_c_mol.GetAtoms() if atom.GetAtomicNum() == 6)
    prod_o_count = sum(1 for atom in product_c_mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if not (prod_c_count == start_c_count + 1 and \
            prod_h_count == start_h_count + 4 and \
            prod_o_count == start_o_count + 1):
        return (f"Incorrect Atom Balance: The product formula {product_c_formula} is not consistent "
                f"with the addition of 'CH4O' to the starting material {start_formula}.")

    # 3b: Check if the less reactive double bond is preserved.
    if not product_c_mol.HasSubstructMatch(conjugated_db_pattern):
        return ("Incorrect Product Structure: The product is missing the α,β-unsaturated ketone "
                "moiety ([CH2]=[CH]-C=O), which was expected to be unreactive.")

    # 3c: Check if the more reactive double bond has been transformed correctly.
    # The reaction is: -CH=C(Me)2  -> -CH(OH)-C(Me)3
    # We check for the resulting structural feature: a hydroxyl group next to a t-butyl group.
    # SMARTS for -CH(OH)-C(C)(C)C
    product_feature_pattern = Chem.MolFromSmarts('[CH1](O)C(C)(C)C')
    if not product_c_mol.HasSubstructMatch(product_feature_pattern):
        return ("Incorrect Product Structure: The product does not contain the expected "
                "-CH(OH)-C(CH3)3 fragment that would result from the proposed reaction mechanism.")

    # --- Step 4: Final Conclusion ---
    # All checks passed. The structure of answer C is a chemically consistent result of the
    # reaction pathway described in the LLM's reasoning.
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)
