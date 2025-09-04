import re

def check_chemistry_answer():
    """
    Checks the correctness of the selected answer for the given organic chemistry question.

    The function verifies two main points:
    1. The final product (B) must be correct based on the Stork enamine alkylation mechanism.
    2. The sequence of reagents (A) must be chemically logical for a multi-step synthesis.
    """

    # The final answer provided by the LLM to be checked.
    llm_final_answer = "<<<C>>>"

    # Define the options as provided in the question.
    options = {
        'A': {
            'reagents': "(i) LDA, DME (ii) CH3CH2I (iii) H3O+",
            'product': "pentan-2-one + N,N-dimethylethanamine"
        },
        'B': {
            'reagents': "(i) LDA (ii) DME, CH3CH2I, H3O+",
            'product': "pentan-2-one + N,N-dimethylethanamine"
        },
        'C': {
            'reagents': "(i) LDA, DME (ii) CH3CH2I (iii) H3O+",
            'product': "heptan-4-one"
        },
        'D': {
            'reagents': "(i) LDA (ii) DME, CH3CH2I, H3O+",
            'product': "heptan-4-one"
        }
    }

    # --- Chemical Analysis ---
    # 1. Determine the correct product.
    # The reaction is a Stork enamine alkylation of a pentan-2-one derivative.
    # The bulky base LDA favors deprotonation at the less substituted alpha-carbon (C1).
    # Alkylation with ethyl iodide adds a -CH2CH3 group to C1.
    # The resulting ketone is CH3CH2-CH2-C(=O)-CH2CH2CH3.
    # The IUPAC name for this product is heptan-4-one.
    correct_product = "heptan-4-one"

    # 2. Determine the correct reagent sequence.
    # The reaction must be stepwise: 1. Base, 2. Alkylating agent, 3. Hydrolysis.
    # A sequence that mixes reagents from different steps is incorrect.
    # The correct format should show three distinct steps, indicated by "(i)", "(ii)", and "(iii)".
    def is_sequence_correct(reagent_string):
        return "(i)" in reagent_string and "(ii)" in reagent_string and "(iii)" in reagent_string

    # --- Verification ---
    try:
        # Extract the letter from the answer string, e.g., "C" from "<<<C>>>"
        selected_key = re.search(r'<<<([A-D])>>>', llm_final_answer).group(1)
    except (AttributeError, IndexError):
        return f"Invalid answer format: {llm_final_answer}. Expected format is <<<A>>>."

    if selected_key not in options:
        return f"Invalid option '{selected_key}' selected. Please choose from A, B, C, D."

    selected_option = options[selected_key]

    # Check if the product in the selected option is correct.
    if selected_option['product'] != correct_product:
        return (f"Incorrect. The final answer '{selected_key}' is wrong because the product is incorrect. "
                f"The expected product is '{correct_product}', but the answer states it is '{selected_option['product']}'.")

    # Check if the reagent sequence in the selected option is correct.
    if not is_sequence_correct(selected_option['reagents']):
        return (f"Incorrect. The final answer '{selected_key}' is wrong because the reagent sequence is not chemically sound. "
                f"The reaction must be performed in three distinct steps. The provided sequence '{selected_option['reagents']}' "
                f"incorrectly groups reagents from different steps.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_chemistry_answer()
print(result)