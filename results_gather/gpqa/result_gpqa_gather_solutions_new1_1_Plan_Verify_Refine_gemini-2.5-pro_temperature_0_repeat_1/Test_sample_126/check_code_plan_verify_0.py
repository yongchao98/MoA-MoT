import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def get_molecular_formula(smiles):
    """Calculates the molecular formula from a SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        return f"Error processing SMILES: {e}"

def get_canonical_smiles(smiles):
    """Generates the canonical SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except Exception as e:
        return f"Error processing SMILES: {e}"

def check_answer():
    """
    Checks the correctness of the answer for the Cope rearrangement of 5-butylnona-2,6-diene.
    """
    # --- Define Reactant and Options ---
    # SMILES strings are generated based on IUPAC names.
    reactant = {
        "name": "5-butylnona-2,6-diene",
        "smiles": "CCC=CC(CCCC)C=CCC"
    }
    
    options = {
        "A": {
            "name": "5-ethylundeca-2,6-diene",
            "smiles": "CCCCCC=CC(CC)C=CCC"
        },
        "B": {
            "name": "5-ethyl-4-methyldeca-2,6-diene",
            "smiles": "CCCC=CC(C)C(CC)C=CC"
        },
        "C": {
            "name": "5-ethyl-4-methyldeca-2,6-diene",
            "smiles": "CCCC=CC(C)C(CC)C=CC"
        },
        "D": {
            "name": "4-ethyl-3-methyldeca-1,5-diene",
            "smiles": "CCCCC=CC(CC)C(C)C=C"
        }
    }
    
    # The provided answer to check
    given_answer_key = "D"
    
    # --- Constraint 1: Isomerism ---
    # A rearrangement reaction must produce an isomer of the reactant.
    reactant_formula = get_molecular_formula(reactant["smiles"])
    answer_formula = get_molecular_formula(options[given_answer_key]["smiles"])
    
    if reactant_formula != answer_formula:
        return f"Incorrect: The chosen answer (Option {given_answer_key}) is not an isomer of the reactant. Reactant formula is {reactant_formula}, but answer formula is {answer_formula}."

    # As a sanity check, let's verify all options.
    # C13H24 is the correct formula for the reactant.
    # Option A is C14H26.
    # Option B/C is C13H24.
    # Option D is C13H24.
    if get_molecular_formula(options["A"]["smiles"]) == reactant_formula:
         return f"Internal Check Error: Option A was expected to be a non-isomer but has formula {reactant_formula}."
    
    # --- Constraint 2: Reaction Mechanism Product ---
    # The product must match the structure predicted by the Cope rearrangement mechanism.
    # The step-by-step derivation shows the product is 4-ethyl-3-methyldeca-1,5-diene.
    # Let's verify that the chosen answer (D) matches this derived structure.
    
    # The name of the product derived from first principles
    derived_product_name = "4-ethyl-3-methyldeca-1,5-diene"
    
    # The name of the chosen answer
    chosen_answer_name = options[given_answer_key]["name"]
    
    if chosen_answer_name != derived_product_name:
        return f"Incorrect: The chosen answer (Option {given_answer_key}: {chosen_answer_name}) does not match the product derived from the Cope rearrangement mechanism ({derived_product_name})."

    # --- Constraint 3: Structural Verification with SMILES ---
    # A more robust check is to compare the canonical SMILES strings.
    # The SMILES for the derived product (4-ethyl-3-methyldeca-1,5-diene) is "CCCCC=CC(CC)C(C)C=C".
    derived_product_smiles = "CCCCC=CC(CC)C(C)C=C"
    
    # Get the canonical SMILES for the chosen answer and the derived product
    canonical_answer_smiles = get_canonical_smiles(options[given_answer_key]["smiles"])
    canonical_derived_smiles = get_canonical_smiles(derived_product_smiles)
    
    if canonical_answer_smiles != canonical_derived_smiles:
        return f"Incorrect: The structure of the chosen answer (Option {given_answer_key}) does not match the structure of the product derived from the Cope rearrangement mechanism. Canonical SMILES do not match: {canonical_answer_smiles} vs {canonical_derived_smiles}."

    # --- Constraint 4: Diene Type ---
    # The mechanism predicts a 1,5-diene product. Let's check the diene type of the answer.
    # 4-ethyl-3-methyldeca-1,5-diene is by name a 1,5-diene.
    # Options A, B, C are 2,6-dienes. This is inconsistent with the mechanism.
    if "1,5-diene" not in chosen_answer_name:
        return f"Incorrect: The Cope rearrangement mechanism predicts a 1,5-diene product, but the chosen answer (Option {given_answer_key}) is a {chosen_answer_name.split('-')[-1]}."

    return "Correct"

# Run the check
try:
    result = check_answer()
    print(result)
except ImportError:
    print("RDKit library not found. Please install it using 'pip install rdkit'. This check cannot be performed without it.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")