import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer based on chemical principles.
    """
    question_options = {
        "A": "(2R,3R,4S)-2-benzyl-4-hydroxy-2-methyl-3-phenylcyclohexan-1-one",
        "B": "(1S,2S,4S)-1-(benzyloxy)-2-methyl-1,2,3,4-tetrahydro-[1,1'-biphenyl]-4-ol",
        "C": "(2S,3S,4S)-2-benzyl-4-hydroxy-2-methyl-3-phenylcyclohexan-1-one",
        "D": "(2S,3R,4S,6S)-2-benzyl-4-hydroxy-6-methyl-3-phenylcyclohexan-1-one"
    }
    
    llm_answer = "D"

    # --- Define Correct Chemical Properties based on reaction analysis ---
    # Constraint 1: The backbone must be a cyclohexanone derivative.
    correct_backbone = "cyclohexan-1-one"
    
    # Constraint 2: The methylation must occur at C6. This is the most critical constraint.
    # The C2 carbon in product 2 is quaternary and has no protons for LDA to remove.
    correct_methylation_site_pattern = "6-methyl"
    incorrect_methylation_site_pattern = "2-methyl"
    
    # Constraint 3: The stereochemistry must be (2S,3R,4S,6S).
    correct_stereochemistry_pattern = "(2S,3R,4S,6S)"

    # --- Analyze the LLM's chosen option ---
    chosen_option_text = question_options.get(llm_answer)

    if not chosen_option_text:
        return f"Invalid answer choice: '{llm_answer}'. Not in options A, B, C, D."

    # Check against constraints
    if correct_backbone not in chosen_option_text:
        return (f"Incorrect. The provided answer '{llm_answer}' is wrong because its molecular backbone is not a cyclohexanone. "
                f"The reaction sequence produces a cyclohexan-1-one derivative, but the answer describes '{chosen_option_text}'.")

    if incorrect_methylation_site_pattern in chosen_option_text:
        return (f"Incorrect. The provided answer '{llm_answer}' is wrong because it places the methyl group at the C2 position. "
                "This is chemically impossible, as the C2 carbon has no protons to be removed by LDA in step 3. "
                "The methylation must occur at the C6 position.")

    if correct_methylation_site_pattern not in chosen_option_text:
        # This case handles options that don't have a methyl group or have it elsewhere.
        return (f"Incorrect. The provided answer '{llm_answer}' is wrong because the methyl group is not at the C6 position, "
                "which is the only possible site for methylation by LDA in step 3.")

    if correct_stereochemistry_pattern not in chosen_option_text:
        # Find the stereochemistry in the chosen option to provide a more detailed error.
        found_stereo = re.search(r'\((.*?)\)', chosen_option_text)
        found_stereo_str = found_stereo.group(0) if found_stereo else "not found"
        return (f"Incorrect. The provided answer '{llm_answer}' has the correct regiochemistry but incorrect stereochemistry. "
                f"Expected {correct_stereochemistry_pattern} but got {found_stereo_str}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)