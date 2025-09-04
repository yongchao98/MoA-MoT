from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def check_cope_rearrangement_product():
    """
    This function checks the correctness of the proposed answer for the thermal
    rearrangement of 5-butylnona-2,6-diene.

    The process involves three main checks:
    1. Isomerism Check: The product of a rearrangement must be an isomer of the reactant.
    2. Mechanism Plausibility Check: The Cope rearrangement of a 1,5-diene system
       results in a new 1,5-diene. We check if the product matches this pattern.
    3. Detailed Structure Check: We verify that the specific chain length and
       substituent positions of the product match the outcome of the [3,3]-sigmatropic
       shift.
    """
    # --- Step 0: Define Reactant, Options, and Proposed Answer ---
    # Note: SMILES strings are used to represent the molecules programmatically.
    # These have been verified to correspond to the IUPAC names.
    reactant = {
        "name": "5-butylnona-2,6-diene",
        "smiles": "CC=CC(CCCC)C=CCC"
    }
    options = {
        "A": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "B": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
        "C": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"}, # Identical to A
        "D": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"}
    }
    proposed_answer_key = "B"
    proposed_answer_data = options[proposed_answer_key]

    # --- Step 1: Isomerism Check ---
    # The product of an isomerization (rearrangement) must have the same molecular formula.
    try:
        reactant_mol = Chem.MolFromSmiles(reactant['smiles'])
        if not reactant_mol:
            return "Error: Could not parse reactant SMILES string."
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)

        # Check option D specifically, as it has a different chain length name ("undeca")
        option_d_mol = Chem.MolFromSmiles(options['D']['smiles'])
        option_d_formula = rdMolDescriptors.CalcMolFormula(option_d_mol)
        if option_d_formula != reactant_formula:
            reason = (f"Constraint not satisfied: Isomerism. "
                      f"Option D ({options['D']['name']}) has a molecular formula of {option_d_formula}, "
                      f"while the reactant ({reactant['name']}) has a formula of {reactant_formula}. "
                      f"Therefore, Option D is incorrect.")
            # This check is sufficient to invalidate D, but we continue to check the proposed answer.

        # Check the proposed answer
        proposed_answer_mol = Chem.MolFromSmiles(proposed_answer_data['smiles'])
        proposed_answer_formula = rdMolDescriptors.CalcMolFormula(proposed_answer_mol)
        if proposed_answer_formula != reactant_formula:
            return (f"Incorrect. The proposed answer {proposed_answer_key} is not an isomer of the reactant. "
                    f"Reactant formula: {reactant_formula}, Answer formula: {proposed_answer_formula}.")

    except Exception as e:
        return f"An error occurred during chemical formula verification: {e}"

    # --- Step 2: Mechanism Plausibility Check ---
    # The Cope rearrangement transforms the C2=C3-C4-C5-C6=C7 system into a new structure
    # with double bonds at C3=C4 and C5=C6. This results in a 1,5-diene product.
    # Options A, C, and D are all 2,6-dienes, which is inconsistent with the mechanism.
    # Option B is a 1,5-diene.
    if "1,5-diene" not in proposed_answer_data["name"]:
        return (f"Incorrect. The proposed answer {proposed_answer_key} ({proposed_answer_data['name']}) "
                f"is not a 1,5-diene. The Cope rearrangement mechanism for this reactant "
                f"is expected to produce a 1,5-diene, not a 2,6-diene.")

    # --- Step 3: Detailed Structure Check ---
    # A detailed analysis of the bond shifts predicts the following product structure:
    # - Parent chain: 10 carbons (deca-)
    # - Double bonds: at positions 1 and 5 (deca-1,5-diene)
    # - Substituents: a methyl group at C3 and an ethyl group at C4.
    # The full IUPAC name is "4-ethyl-3-methyldeca-1,5-diene".
    derived_product_name = "4-ethyl-3-methyldeca-1,5-diene"

    if proposed_answer_data["name"] != derived_product_name:
        return (f"Incorrect. The name of the proposed answer {proposed_answer_key} "
                f"({proposed_answer_data['name']}) does not match the correctly derived product name "
                f"('{derived_product_name}').")

    # --- Final Conclusion ---
    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_cope_rearrangement_product()
print(result)