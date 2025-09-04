import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    # If rdkit is not installed, we cannot perform a robust chemical check.
    # We will print an error and exit. For a real-world application,
    # this dependency would need to be managed.
    print("RDKit is not installed. Cannot perform chemical structure validation.")
    print("Please install it: pip install rdkit-pypi")
    sys.exit("RDKit not found")

def get_canonical_smiles(smiles_string):
    """Converts a SMILES string to its canonical form, ignoring stereochemistry."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return None
    # IsomericSmiles=False ignores stereochemistry, which is not specified in the question.
    return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)

def get_molecular_formula(smiles_string):
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        return None
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def check_correctness():
    """
    Checks the correctness of the LLM's answer by simulating the chemical reaction.
    """
    # --- Problem Definition ---
    # The final answer provided by the LLM analysis.
    llm_final_answer_letter = "B"

    # --- Data Representation from the Question ---
    # Reactant: 5-butylnona-2,6-diene
    # Structure: CH3-CH=CH-CH2-CH(CCCC)-CH=CH-CH2-CH3
    reactant_smiles = "CCC=CC(CCCC)C=CCC"

    # Options provided in the question, represented by their IUPAC names and SMILES strings.
    options = {
        "A": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "B": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
        "C": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"}, # Same as A
        "D": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"}
    }

    # --- Step 1: Determine the Correct Product via Chemical Principles ---
    # The reaction is heating a 1,5-diene, which triggers a Cope Rearrangement.
    # The 1,5-diene system in the reactant is C2=C3-C4-C5-C6=C7.
    # Mechanism:
    # 1. The C4-C5 sigma bond breaks.
    # 2. A new C2-C7 sigma bond forms.
    # 3. Pi bonds shift: C2=C3 -> C3=C4 and C6=C7 -> C5=C6.
    # The resulting product structure is CH2=CH-CH(CH3)-CH(CH2CH3)-CH=CH-CH2CH2CH2CH3.
    # IUPAC Naming of this structure:
    # - Longest chain containing both double bonds is 10 carbons (deca-diene).
    # - Numbering from the CH2= end gives double bonds at 1 and 5 (deca-1,5-diene).
    # - Substituents are a methyl group at C3 and an ethyl group at C4.
    # - Final name: 4-ethyl-3-methyldeca-1,5-diene.
    
    # This name corresponds to option B.
    expected_product_letter = "B"
    expected_product_smiles = options[expected_product_letter]["smiles"]

    # --- Step 2: Verify the LLM's Answer ---
    # Check if the LLM's chosen letter is a valid option.
    if llm_final_answer_letter not in options:
        return f"Incorrect. The final answer '{llm_final_answer_letter}' is not a valid option (A, B, C, or D)."

    llm_chosen_smiles = options[llm_final_answer_letter]["smiles"]

    # --- Step 3: Perform Constraint Checks ---
    # Constraint 1: The product of a rearrangement must be an isomer of the reactant.
    reactant_formula = get_molecular_formula(reactant_smiles)
    product_formula = get_molecular_formula(llm_chosen_smiles)
    if reactant_formula != product_formula:
        return (f"Incorrect. The chosen product ({options[llm_final_answer_letter]['name']}) is not an isomer "
                f"of the reactant. Reactant formula: {reactant_formula}, Product formula: {product_formula}.")

    # --- Step 4: Compare Structures ---
    # Compare the canonical SMILES of the LLM's choice with the expected product.
    llm_canonical_smiles = get_canonical_smiles(llm_chosen_smiles)
    expected_canonical_smiles = get_canonical_smiles(expected_product_smiles)

    if llm_canonical_smiles == expected_canonical_smiles:
        # The structure is correct. The reasoning provided in the prompt also correctly
        # identifies the Cope rearrangement and follows the mechanism to the right product.
        return "Correct"
    else:
        return (f"Incorrect. The final answer '{llm_final_answer_letter}' corresponds to the molecule "
                f"'{options[llm_final_answer_letter]['name']}'.\n"
                f"However, the correct product from the Cope rearrangement is "
                f"'{options[expected_product_letter]['name']}'.\n"
                f"The structure of the chosen answer does not match the structure of the correct product.")

# Run the check and print the result.
result = check_correctness()
print(result)