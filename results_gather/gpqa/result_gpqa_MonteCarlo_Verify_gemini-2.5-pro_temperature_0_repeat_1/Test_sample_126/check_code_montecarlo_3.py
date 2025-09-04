from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_cope_rearrangement_product():
    """
    Checks the correctness of the LLM's answer for the product of heating 5-butylnona-2,6-diene.

    The function will:
    1. Identify the reactant and the expected reaction type (Cope rearrangement).
    2. Determine the correct product based on the [3,3]-sigmatropic shift mechanism.
    3. Compare the derived correct product with the options provided.
    4. Evaluate if the LLM's answer (C) is correct based on chemical principles.
    """

    # --- 1. Define Reactant and Options ---
    # The reactant is 5-butylnona-2,6-diene.
    # The reaction condition "heated" on a 1,5-diene strongly implies a Cope Rearrangement.
    reactant_name = "5-butylnona-2,6-diene"
    reactant_smiles = "CCC=CC(CCCC)C=CCC"

    options = {
        "A": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
        "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCC=CC(C)C(CC)C=CC"},
        "C": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"},
        "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCC=CC(C)C(CC)C=CC"} # Duplicate of B
    }
    llm_answer_key = "C"

    # --- 2. Verify Isomerism (Constraint Check) ---
    # A rearrangement reaction must result in an isomer.
    try:
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        if not reactant_mol: return f"Error: Invalid SMILES for reactant '{reactant_name}'."
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)

        for key, data in options.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            if not mol: return f"Error: Invalid SMILES for option {key}."
            formula = rdMolDescriptors.CalcMolFormula(mol)
            if formula != reactant_formula:
                return (f"Incorrect. A product must be an isomer of the reactant. "
                        f"Reactant formula: {reactant_formula}, Option {key} formula: {formula}.")
    except Exception as e:
        return f"An error occurred during chemical validation: {e}"

    # --- 3. Determine the Correct Product via Mechanism ---
    # The Cope rearrangement is a [3,3]-sigmatropic shift. For the acyclic reactant:
    # Reactant: Me-CH=CH-CH2-CH(Bu)-CH=CH-Et
    # We number the 6 atoms of the pericyclic system:
    # Me-C(1)=C(2)-C(3)H2-C(4)H(Bu)-C(5)=C(6)-Et
    # The reaction breaks the C(3)-C(4) sigma bond and forms a new C(1)-C(6) sigma bond.
    # The pi bonds shift: C(1)=C(2) -> C(2)=C(3) and C(5)=C(6) -> C(4)=C(5).
    # The resulting product structure is: CH2=CH-CH(Me)-CH(Et)-CH=CH-Bu
    # Let's write the full structure: CH2=CH-CH(CH3)-CH(CH2CH3)-CH=CH-CH2CH2CH2CH3
    # The IUPAC name for this product is 4-ethyl-3-methyldeca-1,5-diene.
    # This corresponds to Option A.

    correct_product_name = "4-ethyl-3-methyldeca-1,5-diene"
    correct_option_key = "A"

    # --- 4. Final Verdict ---
    if llm_answer_key == correct_option_key:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer_key}, but the correct product from the reaction is Option {correct_option_key} ({correct_product_name}).\n"
                f"Reasoning: Heating a 1,5-diene like {reactant_name} causes a Cope rearrangement. This is a [3,3]-sigmatropic shift that, in this case, converts one acyclic diene into another.\n"
                f"The product of this specific rearrangement is {correct_product_name}, which is Option A.\n"
                f"The LLM's answer, Option C ({options['C']['name']}), cannot be formed through this standard, expected mechanism. The LLM's reasoning that the Cope product is cyclic is a common but incorrect analysis for this type of acyclic system, and its subsequent choice by elimination is based on this flawed chemical premise.")

# Execute the check and print the result
result = check_cope_rearrangement_product()
print(result)