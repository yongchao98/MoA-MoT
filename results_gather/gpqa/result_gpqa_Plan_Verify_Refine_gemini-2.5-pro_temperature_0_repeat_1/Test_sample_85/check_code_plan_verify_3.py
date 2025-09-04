import sys
from io import StringIO

def check_answer_correctness():
    """
    This function checks the correctness of the provided LLM's answer to a chemistry question.
    It codifies the chemical principles involved in the reactions to verify the logical deduction.
    """
    
    # Store the original stdout to restore it later
    original_stdout = sys.stdout
    # Redirect stdout to a string buffer to capture prints
    captured_output = StringIO()
    sys.stdout = captured_output

    # --- Step 1: Define the chemical principles and problem constraints ---

    # Principle 1: Chemoselectivity of reducing agents
    # LiBH4 selectively reduces esters over carboxylic acids.
    # BH3 selectively reduces carboxylic acids over esters.
    chemoselectivity = {
        'LiBH4': {'reduces': 'ester', 'spares': 'carboxylic_acid'},
        'BH3': {'reduces': 'carboxylic_acid', 'spares': 'ester'}
    }

    # Principle 2: Stereochemistry
    # The reactions occur at the functional groups (C1-acid and C5-ester).
    # The chiral center is at C3. The bonds to the chiral center are not broken or formed.
    # Therefore, the stereochemistry (R/S configuration) is retained.
    stereochemistry_retained = True

    # Problem Statement:
    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    product_stereochemistry = {
        'A': 'R',
        'B': 'S'
    }

    # --- Step 2: Deduce the required starting material stereochemistry based on the principles ---

    # Deduction for Reaction A
    reagent_A = 'LiBH4'
    product_A_config = product_stereochemistry['A']
    # Since stereochemistry is retained, the starting material must have the same configuration as the product.
    required_A_config = product_A_config if stereochemistry_retained else ('S' if product_A_config == 'R' else 'R')

    # Deduction for Reaction B
    reagent_B = 'BH3'
    product_B_config = product_stereochemistry['B']
    # Since stereochemistry is retained, the starting material must have the same configuration as the product.
    required_B_config = product_B_config if stereochemistry_retained else ('S' if product_B_config == 'R' else 'R')

    # --- Step 3: Define the options and the LLM's answer ---
    options = {
        'A': {'A': 'S', 'B': 'R'},
        'B': {'A': 'R', 'B': 'S'},
        'C': {'A': 'S', 'B': 'S'},
        'D': {'A': 'R', 'B': 'R'}
    }
    
    llm_answer_choice = 'B' # The LLM's final answer is <<<B>>>
    llm_answer_details = options.get(llm_answer_choice)

    # --- Step 4: Verify the LLM's answer against the deduced requirements ---
    
    errors = []
    
    # Check the reasoning for starting material A
    if llm_answer_details['A'] != required_A_config:
        error_msg = (f"Constraint violated for starting material A. "
                     f"The reaction with LiBH4 produces an (R)-product. "
                     f"Since the reaction retains stereochemistry, starting material A must be (R). "
                     f"The selected answer claims A is ({llm_answer_details['A']}).")
        errors.append(error_msg)

    # Check the reasoning for starting material B
    if llm_answer_details['B'] != required_B_config:
        error_msg = (f"Constraint violated for starting material B. "
                     f"The reaction with BH3 produces an (S)-product. "
                     f"Since the reaction retains stereochemistry, starting material B must be (S). "
                     f"The selected answer claims B is ({llm_answer_details['B']}).")
        errors.append(error_msg)

    # Restore stdout
    sys.stdout = original_stdout

    # --- Step 5: Return the final result ---
    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Execute the check and print the result
result = check_answer_correctness()
print(result)