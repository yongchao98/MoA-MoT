import sys

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by simulating the chemical reaction.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return ("Skipping check: The 'rdkit' library is required for this verification. "
                "Please install it using 'pip install rdkit'.")

    # --- Step 1: Define the problem based on the provided text ---

    # Reactant from the question
    reactant_name = "5-butylnona-2,6-diene"
    # SMILES representation of the reactant: CH3-CH=CH-CH2-CH(C4H9)-CH=CH-CH2-CH3
    reactant_smiles = "CCC=CC(CCCC)CC=CC"

    # Options from the question
    options = {
        "A": "5-ethyl-4-methyldeca-2,6-diene",
        "B": "4-ethyl-3-methyldeca-1,5-diene",
        "C": "5-ethyl-4-methyldeca-2,6-diene", # Same as A
        "D": "5-ethylundeca-2,6-diene"
    }
    
    # SMILES representations for each option
    options_smiles = {
        "A": "CCCC=CC(CC)C(C)C=CC",
        "B": "CCCCC=CC(CC)C(C)C=C",
        "C": "CCCC=CC(CC)C(C)C=CC",
        "D": "CCCCC=CC(CC)CC=CC"
    }

    # The final answer provided by the LLM to be checked
    llm_final_choice = "B"

    # --- Step 2: Simulate the Cope Rearrangement ---

    # Create an RDKit molecule object for the reactant
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    if not reactant_mol:
        return "Error: Failed to parse the reactant's SMILES string."

    # Define the Cope rearrangement ([3,3]-sigmatropic shift) using a reaction SMARTS string.
    # This rule finds a 1,5-diene and rearranges the bonds accordingly.
    # Pattern: [1]=[2]-[3]-[4]-[5]=[6] >> [2]=[3]-[4]=[5]-[6]-[1]
    rxn_smarts = "[*:1]=[*:2]-[*:3]-[*:4]-[*:5]=[*:6]>>[*:2]=[*:3]-[*:4]=[*:5]-[*:6]-[*:1]"
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    except Exception as e:
        return f"Error: Invalid reaction SMARTS definition. {e}"

    # Run the reaction on the reactant molecule
    products = rxn.RunReactants((reactant_mol,))

    if not products or not products[0]:
        return "Error: The Cope rearrangement simulation failed to produce a product."

    # Get the first product molecule and sanitize it
    product_mol = products[0][0]
    try:
        Chem.SanitizeMol(product_mol)
    except Exception as e:
        return f"Error: The generated product molecule is chemically invalid. {e}"

    # --- Step 3: Compare and Verify the Answer ---

    # To compare molecules, we use their canonical SMILES representation, which is a unique
    # string for a given molecular structure (ignoring stereochemistry).
    actual_product_canonical_smiles = Chem.MolToSmiles(product_mol, canonical=True)

    # Find which option matches the calculated product
    correct_option_letter = None
    for letter, smiles in options_smiles.items():
        option_mol = Chem.MolFromSmiles(smiles)
        if option_mol:
            option_canonical_smiles = Chem.MolToSmiles(option_mol, canonical=True)
            if actual_product_canonical_smiles == option_canonical_smiles:
                correct_option_letter = letter
                break
    
    if correct_option_letter is None:
        return (f"The calculated product does not match any of the given options. "
                f"Calculated product SMILES: {actual_product_canonical_smiles}")

    # Check if the LLM's choice matches the computationally verified correct option
    if llm_final_choice == correct_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The LLM chose option {llm_final_choice} ({options[llm_final_choice]}), "
                f"but the correct product of the Cope rearrangement is option {correct_option_letter} "
                f"({options[correct_option_letter]}).")

# Execute the check and print the result
result = check_chemistry_answer()
print(result)