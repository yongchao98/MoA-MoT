import re

def check_correctness(question, candidate_answers, final_answer):
    """
    Checks the correctness of the final answer for a multi-step organic synthesis problem.
    The function verifies the key reasoning steps, particularly the regiochemistry and stereochemistry.

    Args:
        question (str): The problem description.
        candidate_answers (str): A string containing all candidate answers.
        final_answer (str): The final chosen answer with detailed reasoning.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Define the ground truth based on established chemical principles ---

    # Step 1: Protection of (S)-4-hydroxycyclohex-2-en-1-one with TBSCl.
    # This is straightforward and preserves the (S) stereocenter at C4.

    # Step 2: Tandem conjugate addition (Ph2CuLi) and alkylation (BnBr).
    # - Conjugate addition of Ph- to C3 occurs anti to the bulky C4-OTBS group.
    #   This establishes a (3R, 4S) trans relationship.
    # - Alkylation with BnBr at C2 occurs anti to the newly added C3-Ph group.
    #   This establishes a (2S, 3R) trans relationship.
    # - Expected stereochemistry of Product 2: (2S,3R,4S).

    # Step 3: Methylation with LDA and CH3I.
    # - LDA is a strong, bulky base. At low temperatures, it forms the kinetic enolate.
    # - The kinetic enolate is formed by deprotonating the less sterically hindered alpha-carbon.
    # - In Product 2, C2 is tertiary and heavily hindered. C6 is secondary and less hindered.
    # - Therefore, methylation MUST occur at C6. This is the most critical point for distinguishing options.
    expected_methylation_site = "C6"
    # - The stereochemistry of methylation at C6 is also controlled, leading to a (6S) configuration.
    expected_final_stereochem = "(2S,3R,4S,6S)"

    # Step 4: Deprotection with aqueous HCl.
    # - This simply removes the TBS group, regenerating the alcohol.

    # Final Product:
    expected_final_name = "(2S,3R,4S,6S)-2-benzyl-4-hydroxy-6-methyl-3-phenylcyclohexan-1-one"
    expected_final_option = "D"

    # --- Analyze the provided final answer ---

    # 1. Check the final chosen option letter.
    match = re.search(r'<<<([A-D])>>>', final_answer)
    if not match:
        return "Incorrect: The final answer does not contain a choice in the format <<<X>>>."
    chosen_option = match.group(1)

    # 2. Check the most critical reasoning point: the site of methylation.
    # The analysis must correctly identify that LDA deprotonates at C6.
    methylation_at_c6_reasoning = re.search(r'(deprotonate|methylation|methyl group).* at .*C6', final_answer, re.IGNORECASE)
    methylation_at_c2_reasoning = re.search(r'(deprotonate|methylation|methyl group).* at .*C2', final_answer, re.IGNORECASE)

    if not methylation_at_c6_reasoning:
        if methylation_at_c2_reasoning:
            return ("Incorrect: The analysis incorrectly concludes that methylation occurs at C2. "
                    "With the bulky base LDA at low temperature, the kinetic enolate is formed by deprotonating the less hindered C6 position.")
        return ("Incorrect: The analysis does not clearly justify the site of methylation. "
                "Correctly identifying methylation at C6 (due to kinetic enolate formation) is crucial to eliminate options A and B.")

    # 3. Check if the analysis correctly uses this fact to eliminate other options.
    elimination_reasoning = re.search(r'eliminates options A and B.*methylation at C2', final_answer, re.IGNORECASE)
    if not elimination_reasoning:
        return ("Incorrect: The analysis correctly identifies methylation at C6 but fails to explicitly state that this "
                "is the reason why options A and B (which show methylation at C2) are wrong.")

    # 4. Check if the derived final product name and stereochemistry match the expected result.
    if expected_final_name not in final_answer:
        return f"Incorrect: The derived final product name in the analysis does not match the expected correct name '{expected_final_name}'."

    # 5. Check if the final chosen option is consistent with the correct reasoning.
    if chosen_option != expected_final_option:
        return f"Incorrect: The final chosen option is <<<{chosen_option}>>>, but the correct option based on the chemical principles is <<<{expected_final_option}>>>."

    # If all checks pass, the answer is correct.
    return "Correct"
