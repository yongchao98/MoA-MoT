import sys
from io import StringIO

# It's recommended to install rdkit: pip install rdkit
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError:
    print("RDKit not found. Please install it using 'pip install rdkit'")
    # Create dummy functions if rdkit is not available to avoid crashing the script
    # The script will fail the check but will not raise an ImportError
    class Chem:
        @staticmethod
        def MolFromSmiles(s): return None
        @staticmethod
        def MolToSmiles(m, isomericSmiles=True, canonical=True): return None
    class rdMolDescriptors:
        @staticmethod
        def CalcMolFormula(m): return None


def get_canonical_smiles(smiles: str) -> str:
    """Converts a SMILES string to its canonical form."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    # The isomericSmiles=True is important to distinguish stereoisomers if any were present
    return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

def get_molecular_formula(smiles: str) -> str:
    """Calculates the molecular formula from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return rdMolDescriptors.CalcMolFormula(mol)

def check_correctness():
    """
    Checks the correctness of the final answer by simulating the chemical reasoning.
    1. Defines the reactant and options using their chemical structures (SMILES).
    2. Determines the expected product based on the Cope rearrangement mechanism.
    3. Verifies that the expected product is an isomer of the reactant.
    4. Verifies that the expected product matches the expected product type (a 1,5-diene).
    5. Finds which option corresponds to the correct product structure.
    6. Compares this correct option with the provided final answer.
    7. Performs sanity checks on the incorrect options.
    """
    try:
        # Step 1: Define the problem data based on the question and the final answer's analysis.
        # Using hardcoded SMILES strings for robustness, avoiding reliance on external name-to-SMILES services.
        reactant_name = "5-butylnona-2,6-diene"
        reactant_smiles = "CC=CCC(CCCC)C=CCC"

        # The options as listed in the final answer's analysis section
        options = {
            "A": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCC=CC(C)C(CC)C=CCC"},
            "B": {"name": "4-ethyl-3-methyldeca-1,5-diene", "smiles": "C=CC(C)C(CC)C=CCCCC"},
            "C": {"name": "5-ethylundeca-2,6-diene", "smiles": "CCCC=CC(CC)C=CCCC"},
            "D": {"name": "5-ethyl-4-methyldeca-2,6-diene", "smiles": "CCC=CC(C)C(CC)C=CCC"} # Same as A
        }
        
        # The final answer provided by the LLM to be checked
        llm_final_answer = "B"

        # Step 2: Determine the expected product based on chemical principles.
        # The Cope rearrangement of 5-butylnona-2,6-diene yields 4-ethyl-3-methyldeca-1,5-diene.
        expected_product_name = "4-ethyl-3-methyldeca-1,5-diene"
        expected_product_smiles = "C=CC(C)C(CC)C=CCCCC"

        # Step 3: Verify key constraints.
        
        # Constraint 1: Isomerism. The product of a rearrangement must be an isomer of the reactant.
        reactant_formula = get_molecular_formula(reactant_smiles)
        expected_product_formula = get_molecular_formula(expected_product_smiles)
        if not reactant_formula or not expected_product_formula:
             return "Could not generate molecular formulas. RDKit might be missing or SMILES are invalid."
        if reactant_formula != expected_product_formula:
            return f"Reason: The chemically correct product ({expected_product_name}, formula {expected_product_formula}) is not an isomer of the reactant ({reactant_name}, formula {reactant_formula})."

        # Constraint 2: Reaction product type. The Cope rearrangement should yield a 1,5-diene.
        if "1,5-diene" not in expected_product_name:
            return f"Reason: The chemically correct product '{expected_product_name}' is not a 1,5-diene, which is the expected product type for a Cope rearrangement."

        # Step 4: Match the derived product to the options by comparing canonical SMILES.
        expected_canon_smiles = get_canonical_smiles(expected_product_smiles)
        
        correct_option_letter = None
        for letter, data in options.items():
            option_canon_smiles = get_canonical_smiles(data["smiles"])
            if option_canon_smiles == expected_canon_smiles:
                correct_option_letter = letter
                # We can break here, but iterating through all ensures no duplicates for the correct answer
        
        if correct_option_letter is None:
            return f"Reason: The chemically correct product '{expected_product_name}' was not found among the provided options."

        # Step 5: Compare the identified correct option with the LLM's final answer.
        if correct_option_letter != llm_final_answer:
            return f"Reason: The correct product is '{expected_product_name}', which corresponds to option {correct_option_letter}, but the provided answer was {llm_final_answer}."

        # Step 6: Perform sanity checks on the other options to confirm the analysis.
        # Option C should not be an isomer.
        option_c_formula = get_molecular_formula(options["C"]["smiles"])
        if option_c_formula == reactant_formula:
            return f"Reason: The analysis is flawed. Option C ({options['C']['name']}) was expected to be a non-isomer, but its formula ({option_c_formula}) matches the reactant's."
        
        # Options A/D are isomers but are 2,6-dienes, not the 1,5-diene product.
        if "2,6-diene" not in options["A"]["name"]:
            return f"Reason: The analysis is flawed. Option A ({options['A']['name']}) was expected to be a 2,6-diene, but its name does not reflect this."

        # If all checks pass, the LLM's reasoning and final answer are correct.
        return "Correct"

    except Exception as e:
        # Catch any unexpected errors during execution.
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)