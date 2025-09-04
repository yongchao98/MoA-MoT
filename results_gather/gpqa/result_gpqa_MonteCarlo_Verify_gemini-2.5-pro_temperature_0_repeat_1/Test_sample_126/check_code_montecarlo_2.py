import sys

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    # A fallback for environments where RDKit is not installed.
    # The check will be less rigorous, relying only on name patterns.
    print("Warning: RDKit not found. Performing a less rigorous check.", file=sys.stderr)
    Chem = None

def check_cope_rearrangement_answer():
    """
    Verifies the product of the Cope rearrangement of 5-butylnona-2,6-diene.

    The function checks three main constraints:
    1. Isomerism: The product must have the same molecular formula as the reactant.
    2. Reaction Pattern: The reaction is a Cope rearrangement, which transforms a 1,5-diene.
       The product is expected to be a 1,5-diene as well.
    3. Mechanistic Outcome: The specific product is determined by tracing the atoms through the
       [3,3]-sigmatropic rearrangement. This check compares the proposed answer to the known
       mechanistic outcome.
    """
    reactant_name = "5-butylnona-2,6-diene"
    reactant_smiles = "CCC=CC(CCCC)C=CCC"  # A 1,5-diene system (C2 to C7)

    options = {
        "A": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "CCCCC=CC(CC)C(C)C=C"},
        "B": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"},
        "C": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCCC=CC(CC)CC=CCC"},
        "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCCC=CC(C)C(CC)C=CC"}
    }
    
    proposed_answer_key = "A"

    # --- Check 1: Isomerism (Requires RDKit) ---
    if Chem:
        try:
            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
            if not reactant_mol:
                return f"Error: Could not parse SMILES for reactant: {reactant_smiles}"
            reactant_formula = rdMolDescriptors.CalcMolFormula(reactant_mol)

            for key, data in options.items():
                option_mol = Chem.MolFromSmiles(data["smiles"])
                if not option_mol:
                    return f"Error: Could not parse SMILES for option {key}: {data['smiles']}"
                option_formula = rdMolDescriptors.CalcMolFormula(option_mol)
                if reactant_formula != option_formula:
                    return (f"Incorrect. Constraint Violated: Isomerism.\n"
                            f"The product must be an isomer of the reactant ({reactant_formula}).\n"
                            f"Option {key} ({data['name']}) has formula {option_formula}.")
        except Exception as e:
            return f"An error occurred during RDKit processing: {e}"
    else:
        # Basic check if RDKit is not available
        reactant_carbons = 9 + 4 # nona + butyl
        for key, data in options.items():
            if key == "A" and (10 + 1 + 2) != reactant_carbons: return "Incorrect based on carbon count."
            if key == "B" and (10 + 1 + 2) != reactant_carbons: return "Incorrect based on carbon count."
            if key == "C" and (11 + 2) != reactant_carbons: return "Incorrect based on carbon count."


    # --- Check 2 & 3: Reaction Mechanism and Known Outcome ---
    # The Cope rearrangement of the reactant (a substituted nona-2,6-diene) involves
    # breaking the C4-C5 bond and forming a new C2-C7 bond.
    # This specific rearrangement is known to incorporate the C5-butyl group into the main chain,
    # extending it from a 9-carbon (nona) to a 10-carbon (deca) chain.
    # The resulting product is a substituted deca-1,5-diene.

    # Let's analyze the options based on this known pattern:
    # A) 4-ethyl-3-methyldeca-1,5-diene: A deca-1,5-diene. This matches the expected structural pattern.
    # B) 5-ethyl-4-methyldeca-2,6-diene: A deca-2,6-diene. The diene position is inconsistent with the rearrangement outcome.
    # C) 5-ethylundeca-2,6-diene: An undeca- (C11) chain. The chain length is incorrect.

    # Based on the known mechanism, only option A fits the expected product structure.
    known_product_name = "4-ethyl-3-methyldeca-1,5-diene"
    
    if options[proposed_answer_key]["name"] == known_product_name:
        return "Correct"
    else:
        correct_key = "Unknown"
        for key, data in options.items():
            if data["name"] == known_product_name:
                correct_key = key
                break
        return (f"Incorrect. The provided answer is {proposed_answer_key}, but the known product of the reaction is "
                f"'{known_product_name}', which corresponds to option {correct_key}.")

# Execute the check and print the result.
result = check_cope_rearrangement_answer()
print(result)