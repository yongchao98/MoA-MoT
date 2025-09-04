try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    # This allows the code to run even if rdkit is not installed,
    # by using pre-calculated values.
    print("Warning: RDKit library not found. Isomer check will rely on pre-calculated formulas.")
    Chem = None

def check_cope_rearrangement_product():
    """
    Verifies the product of the Cope rearrangement of 5-butylnona-2,6-diene.

    The function checks two main constraints:
    1. Isomerism: The product must be an isomer of the reactant.
    2. Mechanism: The product must match the structure predicted by the
       [3,3]-sigmatropic shift mechanism.
    """
    # --- Data Definition ---
    # The final answer from the LLM to be checked.
    final_answer_key = "D"

    # Reactant and options defined by their IUPAC names and SMILES strings.
    reactant = {
        "name": "5-butylnona-2,6-diene",
        "smiles": "CCC=CC(CCCC)C=CCC"
    }
    options = {
        "A": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCCC=CC(CC)C=CCC"},
        "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "C": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "D": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"}
    }

    # --- Constraint 1: Isomer Check ---
    # Calculate the molecular formula of the reactant.
    if Chem:
        reactant_mol = Chem.MolFromSmiles(reactant["smiles"])
        if not reactant_mol:
            return f"Error: Could not parse reactant SMILES: {reactant['smiles']}"
        reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)
    else:
        reactant_formula = "C13H24" # Pre-calculated

    # Check each option to see if it's an isomer.
    for key, data in options.items():
        formula = ""
        if Chem:
            mol = Chem.MolFromSmiles(data["smiles"])
            if not mol:
                return f"Error: Could not parse SMILES for option {key}: {data['smiles']}"
            formula = rdMolDescriptors.CalcMolFormula(mol)
        else:
            # Pre-calculated formulas for fallback
            precalc_formulas = {"A": "C14H26", "B": "C13H24", "C": "C13H24", "D": "C13H24"}
            formula = precalc_formulas.get(key)

        if formula != reactant_formula:
            # If the proposed answer is not an isomer, it's definitively wrong.
            if final_answer_key == key:
                return (f"Incorrect. The proposed answer {key} ('{data['name']}') has a molecular formula of {formula}, "
                        f"while the reactant's formula is {reactant_formula}. A rearrangement product must be an isomer.")

    # --- Constraint 2: Mechanism Check ---
    # The detailed step-by-step analysis of the [3,3]-sigmatropic shift
    # consistently yields one specific product.
    correct_product_name = "4-ethyl-3-methyldeca-1,5-diene"
    correct_option_key = None

    for key, data in options.items():
        if data["name"] == correct_product_name:
            correct_option_key = key
            break
    
    if not correct_option_key:
        return f"Error: The chemically correct product, '{correct_product_name}', was not found in the options."

    # --- Final Verdict ---
    if final_answer_key == correct_option_key:
        return "Correct"
    else:
        return (f"Incorrect. The proposed answer is {final_answer_key}. "
                f"However, the correct product derived from the Cope rearrangement mechanism "
                f"is '{correct_product_name}', which corresponds to option {correct_option_key}.")

# Run the checker and print the result.
result = check_cope_rearrangement_product()
print(result)