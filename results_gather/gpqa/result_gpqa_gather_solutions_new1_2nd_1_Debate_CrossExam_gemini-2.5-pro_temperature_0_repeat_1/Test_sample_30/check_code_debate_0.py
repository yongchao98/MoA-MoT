import re

def check_correctness_of_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer by verifying the chemical reasoning step-by-step.

    Args:
        final_answer_text: The text of the final answer to be checked.

    Returns:
        "Correct" if the answer is correct, or a string explaining the error.
    """
    # --- Ground Truth Data ---
    # The options provided in the original question prompt
    options = {"A": "cs", "B": "c3", "C": "c2h", "D": "d2h"}

    # Define the correct intermediates and products based on sound chemical principles
    product_1_name = "p-nitrotoluene"
    product_2_name = "p-nitrobenzaldehyde"
    
    # The two plausible outcomes for step 3
    product_3_single = {
        "name": "(E)-4-(4-nitrophenyl)but-3-en-2-one",
        "symmetry": "Cs"
    }
    product_3_double = {
        "name": "(1E,4E)-1,5-bis(4-nitrophenyl)penta-1,4-dien-3-one",
        "symmetry": "C2h"
    }

    # --- Analysis of the Provided Answer ---

    # 1. Extract the final letter choice (e.g., <<<C>>>)
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Incorrect. The final answer is not in the required format '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>'."
    
    final_choice_key = match.group(1)
    final_choice_symmetry = options.get(final_choice_key)

    # 2. Verify the logical flow of the chemical analysis in the text
    
    # Check Step 1: Correct identification of Product 1
    if product_1_name.lower() not in final_answer_text.lower():
        return f"Incorrect. The analysis fails to identify Product 1 as {product_1_name}."

    # Check Step 2: Correct identification of Product 2 using "look-ahead" logic
    if product_2_name.lower() not in final_answer_text.lower():
        return f"Incorrect. The analysis fails to identify Product 2 as {product_2_name}."
    if "look-ahead" not in final_answer_text.lower() and "next step" not in final_answer_text.lower():
        return f"Incorrect. The reasoning for choosing {product_2_name} is weak. A strong answer must use 'look-ahead' logic based on the requirements of the subsequent Claisen-Schmidt condensation."

    # Check Step 3: Correct identification of the most plausible Product 3
    # For a symmetry problem, the double condensation product is the most likely intended answer.
    if product_3_double["name"].lower() not in final_answer_text.lower():
        return f"Incorrect. The analysis fails to identify the most plausible final product for a symmetry problem, which is the double condensation product: '{product_3_double['name']}'."

    # Check Step 4: Correct symmetry analysis of the chosen product
    # The text must correctly identify the point group as C2h for the double condensation product.
    # Using regex to find 'C2h' or 'C₂h' followed by 'point group'
    if not re.search(r'c[₂2]h\s+point\s+group', final_answer_text, re.IGNORECASE):
        return f"Incorrect. The symmetry analysis is flawed. It fails to conclude that the double condensation product belongs to the {product_3_double['symmetry']} point group."

    # A high-quality answer should also correctly refute the common D2h mistake.
    if not re.search(r'rules\s+out\s+the\s+d[₂2]h', final_answer_text, re.IGNORECASE):
        return f"Incorrect. The symmetry analysis is incomplete. A robust answer should not only identify the correct point group ({product_3_double['symmetry']}) but also explain why the common incorrect alternative (D2h) is wrong."

    # 3. Final Consistency Check
    # The derived symmetry in the text must match the final letter choice.
    if final_choice_symmetry.lower() != product_3_double["symmetry"].lower():
        return f"Incorrect. There is a contradiction. The analysis correctly derives the point group as {product_3_double['symmetry']}, but the final answer choice is '<<<{final_choice_key}>>>' which corresponds to {final_choice_symmetry}."

    return "Correct"
