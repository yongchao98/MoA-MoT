import re

def check_correctness_of_chemistry_answer(llm_answer_text: str):
    """
    This function checks the correctness of a given LLM's answer for a specific
    organic chemistry question involving the ring-opening of an epoxide.

    It validates the answer based on two key principles:
    1.  Regioselectivity: The nucleophile attacks the less sterically hindered carbon.
    2.  Stereoselectivity: The Sₙ2 attack causes an inversion of configuration, but the
        final R/S descriptor must be determined by re-evaluating CIP priorities.

    Args:
        llm_answer_text: The full text of the LLM's answer, including its reasoning
                         and final choice in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a reason for why the answer is incorrect.
    """

    # --- Ground Truth based on chemical principles and the question's options ---
    # The correct product is (1R,2S,4R,5R)-1,2,4,5-tetramethylcyclohexan-1-ol.
    # This corresponds to option C in the provided list.
    correct_option = "C"
    correct_product_skeleton = "1,2,4,5-tetramethylcyclohexan-1-ol"
    correct_c2_config_str = "(2S)"
    incorrect_c2_config_str = "(2R)"
    incorrect_product_skeleton = "2,2,4,5-tetramethylcyclohexan-1-ol"

    # --- Parse the LLM's Answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer does not provide a final choice in the required format <<<X>>>."
    
    llm_choice = match.group(1)
    # Use the original case for searching for specific terms like (2S)
    llm_reasoning_cased = llm_answer_text
    # Use lowercase for general keyword searches
    llm_reasoning_lower = llm_answer_text.lower()

    # --- Verification Logic ---

    # 1. Check if the final choice is correct. If not, provide a specific reason.
    if llm_choice != correct_option:
        # Reason 1: Incorrect regioselectivity (attack at C1).
        if incorrect_product_skeleton.replace('-', '') in llm_reasoning_lower.replace('-', '') or \
           llm_choice in ['A', 'B']:
            return (f"Incorrect. The final choice is {llm_choice}, but the correct answer is {correct_option}. "
                    f"The answer is based on incorrect regioselectivity. The nucleophilic attack "
                    f"should occur at the less hindered C6, leading to a '{correct_product_skeleton}' skeleton, "
                    f"not a '{incorrect_product_skeleton}' skeleton which would result from an attack at C1.")

        # Reason 2: Incorrect stereoselectivity (the common S -> R error).
        if incorrect_c2_config_str in llm_reasoning_cased or \
           "s inverts to r" in llm_reasoning_lower or \
           llm_choice == 'D':
            return (f"Incorrect. The final choice is {llm_choice}, but the correct answer is {correct_option}. "
                    f"The answer makes a common error in stereochemistry. While the reaction causes geometric "
                    f"inversion at C6, simply flipping the stereodescriptor from (S) to (R) is incorrect due to "
                    f"changes in Cahn-Ingold-Prelog (CIP) priorities. A rigorous 3D analysis shows the product "
                    f"has a {correct_c2_config_str} configuration.")
        
        # Fallback for other errors.
        return f"Incorrect. The final choice is {llm_choice}, but the correct answer is {correct_option}."

    # 2. The final choice is correct (C). Now, verify the reasoning is sound and not contradictory.
    
    # Check for sound regioselectivity reasoning.
    if "c6" not in llm_reasoning_lower or "less hindered" not in llm_reasoning_lower:
        return (f"Incorrect. Although the final choice '{llm_choice}' is correct, the reasoning is flawed. "
                f"It fails to explicitly state that the attack occurs at the less hindered carbon, C6.")

    if correct_product_skeleton.replace('-', '') not in llm_reasoning_lower.replace('-', ''):
        return (f"Incorrect. Although the final choice '{llm_choice}' is correct, the reasoning is flawed. "
                f"It fails to deduce the correct '{correct_product_skeleton}' product skeleton.")

    # Check for sound stereoselectivity reasoning.
    if "inversion" not in llm_reasoning_lower:
        return (f"Incorrect. Although the final choice '{llm_choice}' is correct, the reasoning is flawed. "
                f"It fails to mention the key stereochemical principle of Sₙ2 inversion.")

    # Check for contradiction: choosing C but reasoning for D.
    if incorrect_c2_config_str in llm_reasoning_cased or "s inverts to r" in llm_reasoning_lower:
         return (f"Incorrect. The final choice '{llm_choice}' is correct, but the reasoning is contradictory. "
                 f"It incorrectly states that the configuration at C2 becomes (R) but then selects option C, which has a (2S) configuration.")

    # Check that the correct stereochemistry is actually derived.
    if correct_c2_config_str not in llm_reasoning_cased:
         return (f"Incorrect. Although the final choice '{llm_choice}' is correct, the reasoning is incomplete or flawed. "
                 f"It fails to correctly determine and state that the final stereochemistry at C2 is {correct_c2_config_str}.")

    return "Correct"
