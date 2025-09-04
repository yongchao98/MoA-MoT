import sys
from io import StringIO

def check_answer():
    """
    Checks the correctness of the LLM's answer by simulating the Cope rearrangement
    of 5-butylnona-2,6-diene using RDKit and comparing the result to the given options.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return "Skipping check: RDKit is not installed. Please install it using 'pip install rdkit'."

    # --- Step 1: Define the starting material and the options ---

    # Starting material: 5-butylnona-2,6-diene
    # Structure: CH3-CH=CH-CH2-CH(CCCC)-CH=CH-CH2-CH3
    start_smiles = "CCC=CC(CCCC)CC=CC"
    mol_start = Chem.MolFromSmiles(start_smiles)

    if mol_start is None:
        return "Error: Could not parse the SMILES for the starting material '5-butylnona-2,6-diene'."

    # Options provided in the question
    options = {
        "A": "5-ethyl-4-methyldeca-2,6-diene",
        "B": "5-ethylundeca-2,6-diene",
        "C": "5-ethyl-4-methyldeca-2,6-diene", # Same as A
        "D": "4-ethyl-3-methyldeca-1,5-diene"
    }
    
    # The LLM's final answer
    llm_answer_key = "A"

    # --- Step 2: Define the Cope rearrangement reaction ---

    # A Cope rearrangement is a [3,3]-sigmatropic shift.
    # The pattern C1=C2-C3-C4-C5=C6 rearranges to C2=C3-C4=C5-C6-C1.
    # The C3-C4 bond breaks, and a new C1-C6 bond forms.
    reaction_smarts = "[C:1]=[C:2]-[C:3]-[C:4]-[C:5]=[C:6]>>[C:2]=[C:3]-[C:4]=[C:5]-[C:6]-[C:1]"
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)

    # --- Step 3: Run the reaction ---
    
    products = rxn.RunReactants((mol_start,))

    if not products or not products[0]:
        return "Reaction failed: The Cope rearrangement SMARTS pattern did not match the starting material."

    # We expect one product from this intramolecular rearrangement
    product_mol = products[0][0]
    try:
        Chem.SanitizeMol(product_mol)
    except Exception as e:
        return f"Error in product molecule structure: {e}"

    # --- Step 4: Identify the correct product by comparing SMILES ---

    # Generate a canonical SMILES for the calculated product (ignoring stereochemistry)
    product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=False, canonical=True)

    # Create canonical SMILES for each option for a robust comparison
    option_smiles = {
        "A": Chem.MolToSmiles(Chem.MolFromSmiles("CCCC=CC(CC)C(C)C=CC"), isomericSmiles=False, canonical=True),
        "B": Chem.MolToSmiles(Chem.MolFromSmiles("CCCCC=CC(CC)CC=CC"), isomericSmiles=False, canonical=True),
        "C": Chem.MolToSmiles(Chem.MolFromSmiles("CCCC=CC(CC)C(C)C=CC"), isomericSmiles=False, canonical=True),
        "D": Chem.MolToSmiles(Chem.MolFromSmiles("CCCCC=CC(CC)C(C)C=C"), isomericSmiles=False, canonical=True),
    }

    # Find which option matches the calculated product
    correct_key = None
    for key, smiles in option_smiles.items():
        if product_smiles == smiles:
            correct_key = key
            break
    
    if correct_key is None:
        return (f"The calculated product does not match any of the options.\n"
                f"Calculated product SMILES: {product_smiles}\n"
                f"Options SMILES: {option_smiles}")

    # --- Step 5: Check if the LLM's answer is correct ---

    if llm_answer_key == correct_key:
        return "Correct"
    else:
        return (f"Incorrect. The LLM's answer is '{options[llm_answer_key]}' (Option {llm_answer_key}), "
                f"but the correct product of the Cope rearrangement is '{options[correct_key]}' (Option {correct_key}).\n"
                f"The initial derivation in the LLM's reasoning, which pointed to option D, was correct. "
                f"The final decision to choose option A was based on a flawed premise.")

# Execute the check and print the result
result = check_answer()
print(result)