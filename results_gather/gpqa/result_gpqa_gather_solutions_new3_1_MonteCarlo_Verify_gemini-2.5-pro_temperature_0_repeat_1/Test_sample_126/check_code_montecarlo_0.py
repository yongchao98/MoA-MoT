from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer for the Cope rearrangement question.
    It verifies:
    1. The final chosen option is the correct product.
    2. The reasoning provided is factually accurate.
    """
    # --- Data Setup ---
    # Define molecules from IUPAC names using SMILES strings. These have been manually derived and verified.
    molecules = {
        "Reactant": {"name": "5-butylnona-2,6-diene", "smiles": "CCC=CC(CCCC)C=CCC"},
        "A": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCC=CC(CC)CCC=CC"},
        "B": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
        "C": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(CC)C(C)C=CC"},
    }
    llm_final_choice = "B"

    # --- Analysis ---
    # 1. Determine the correct product based on the reaction mechanism.
    # The reaction is a Cope rearrangement. A step-by-step derivation of the [3,3]-sigmatropic shift
    # on the C2-C7 hexa-1,5-diene system of the reactant yields a single product.
    # The IUPAC name for this product is "4-ethyl-3-methyldeca-1,5-diene".
    correct_product_name = "4-ethyl-3-methyldeca-1,5-diene"
    correct_option_key = "B"

    # 2. Check if the LLM's final choice matches the correct product.
    is_choice_correct = (llm_final_choice == correct_option_key)

    # 3. Check the LLM's reasoning for the stated factual error.
    # The reasoning claims option A is not an isomer and has formula C14H26.
    reasoning_is_flawed = False
    flaw_description = ""
    try:
        reactant_mol = Chem.MolFromSmiles(molecules["Reactant"]["smiles"])
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)  # Expected: C13H24

        mol_A = Chem.MolFromSmiles(molecules["A"]["smiles"])
        formula_A = rdMolDescriptors.CalcMolFormula(mol_A)  # Expected: C13H24

        # If option A is an isomer, the LLM's reasoning is factually incorrect.
        if formula_A == reactant_formula:
            reasoning_is_flawed = True
            flaw_description = (
                f"The reasoning incorrectly states that option A ({molecules['A']['name']}) is not an isomer "
                f"and has a formula of C14H26. The actual formula of option A is {formula_A}, which is an isomer "
                f"of the reactant ({reactant_formula}). The correct reason to eliminate option A is that its structure "
                "is not the one predicted by the Cope rearrangement mechanism."
            )
    except Exception as e:
        return f"An error occurred during molecule processing with RDKit: {e}"

    # --- Final Verdict ---
    if is_choice_correct:
        if reasoning_is_flawed:
            # The final choice is correct, but the justification has a factual error.
            return f"The final choice <<<B>>> is correct, however the reasoning provided is flawed. {flaw_description}"
        else:
            # Both the choice and reasoning are fully correct.
            return "Correct"
    else:
        # The final choice is incorrect.
        return (f"Incorrect. The chosen answer is {llm_final_choice}, but the correct product of the Cope rearrangement "
                f"is '{correct_product_name}', which corresponds to option {correct_option_key}.")

# Execute the check
result = check_chemistry_answer()
print(result)