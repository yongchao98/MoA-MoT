def check_cope_rearrangement_answer():
    """
    Checks the correctness of the answer to the chemistry question by simulating
    the Cope rearrangement using the RDKit library.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return ("RDKit library not found. This check requires RDKit to perform "
                "chemical reaction simulations. Please install it, for example, "
                "using 'pip install rdkit-pypi'.")

    # --- Step 1: Define Reactant and Options ---
    # The molecules are defined by their IUPAC names and corresponding SMILES strings.
    reactant_name = "5-butylnona-2,6-diene"
    reactant_smiles = "CCC=CC(CCCC)C=CCC"

    options = {
        "A": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCC=CC(C)C(CC)C=CCC"},
        "B": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCCC=CC(CC)(C)C=C"},
        "C": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCC=CC(C)C(CC)C=CCC"}, # Same as A
        "D": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"}
    }
    
    # The answer provided by the other LLM.
    llm_answer_key = "B"

    # --- Step 2: Isomer Sanity Check ---
    # A rearrangement must produce an isomer (same molecular formula).
    try:
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if not reactant_mol:
            return f"Error: Could not parse reactant SMILES: {reactant_smiles}"
        reactant_formula = AllChem.CalcMolFormula(reactant_mol)

        answer_mol = Chem.MolFromSmiles(options[llm_answer_key]["smiles"])
        if not answer_mol:
            return f"Error: Could not parse answer SMILES for option {llm_answer_key}: {options[llm_answer_key]['smiles']}"
        answer_formula = AllChem.CalcMolFormula(answer_mol)

        if reactant_formula != answer_formula:
            return (f"Incorrect. The product must be an isomer of the reactant. "
                    f"Reactant formula is {reactant_formula}, but option {llm_answer_key} "
                    f"has formula {answer_formula}.")
    except Exception as e:
        return f"An error occurred during the isomer check with RDKit: {e}"

    # --- Step 3: Define and Run the Cope Rearrangement ---
    # The reaction is defined using SMARTS, a pattern-matching language for molecules.
    # This SMARTS finds a 1,5-diene and performs the [3,3]-sigmatropic shift.
    rxn_smarts = "[C:1]=[C:2]-[C:3]-[C:4]-[C:5]=[C:6]>>[C:2]=[C:3]-[C:4]=[C:5]-[C:6]-[C:1]"
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smarts)
        products = rxn.RunReactants((reactant_mol,))
    except Exception as e:
        return f"An error occurred while running the reaction simulation: {e}"

    if not products:
        return "Incorrect. The reaction simulation did not yield any products for the given reactant."

    # --- Step 4: Verify the Product ---
    # We get a set of possible products. For this reactant, there should be only one.
    # We use canonical SMILES for a standardized, unambiguous comparison.
    product_mols = products[0]
    if not product_mols:
        return "Incorrect. The reaction simulation yielded an empty product set."
        
    predicted_product_mol = product_mols[0]
    Chem.SanitizeMol(predicted_product_mol)
    predicted_canonical_smiles = Chem.MolToSmiles(predicted_product_mol, canonical=True)
    
    # Get the canonical SMILES of the proposed answer molecule.
    answer_canonical_smiles = Chem.MolToSmiles(answer_mol, canonical=True)

    # Compare the predicted product with the proposed answer.
    if predicted_canonical_smiles == answer_canonical_smiles:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer (Option {llm_answer_key}) is not the "
                f"product of a Cope rearrangement.\n"
                f"Reactant: {reactant_name} ({reactant_smiles})\n"
                f"Proposed Answer: {options[llm_answer_key]['name']} ({answer_canonical_smiles})\n"
                f"Predicted Product SMILES: {predicted_canonical_smiles}")

# Execute the check and print the result.
result = check_cope_rearrangement_answer()
print(result)