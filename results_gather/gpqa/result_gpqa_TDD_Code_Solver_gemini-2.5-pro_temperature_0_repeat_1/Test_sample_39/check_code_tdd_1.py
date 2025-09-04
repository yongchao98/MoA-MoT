def check_chemist_slang_answer():
    '''
    This function checks the correctness of the answer to the chemistry slang question.
    It uses a knowledge base of chemistry jargon to determine the correct option and compares it
    to the provided answer.
    '''
    # The key phrase from the question's scenario
    question_phrase = "my compounds are on top of each other"
    
    # The answer provided by the LLM, extracted from the response.
    llm_answer = "A"

    # Define the options and their underlying scientific principles
    options_map = {
        "A": "similar polarities",
        "B": "similar optical rotations",
        "C": "non-covalent bonding",
        "D": "similar boiling points"
    }

    # --- Expert Knowledge Base ---
    # A dictionary mapping common chemist slang to the separation problem it describes.
    # The most common separation technique in a synthetic organic lab is chromatography.
    slang_to_principle = {
        "on top of each other": "similar polarities",  # Refers to chromatography (TLC/column)
        "co-eluting": "similar polarities",  # Refers to column chromatography
        "spots are running together": "similar polarities",  # Refers to TLC (chromatography)
        "can't separate by distillation": "similar boiling points"  # Refers to distillation
    }

    # --- Logic to determine the correct answer ---
    # Find the principle associated with the slang in the question
    # We use .get() to avoid an error if the phrase is not found, though it is present.
    correct_principle = slang_to_principle.get(question_phrase, None)

    if correct_principle is None:
        return f"Constraint not satisfied: The checker's knowledge base does not contain the key phrase '{question_phrase}'."

    # Find which option letter corresponds to the correct principle
    correct_option_letter = None
    for letter, principle in options_map.items():
        if principle == correct_principle:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        return f"Constraint not satisfied: No option corresponds to the identified principle '{correct_principle}'."

    # --- Final Check ---
    # Compare the derived correct answer with the LLM's answer.
    if llm_answer == correct_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The phrase '{question_phrase}' is standard jargon in a chemistry lab for a failed "
                f"chromatographic separation. Chromatography separates compounds based on polarity. Therefore, the "
                f"compounds must have '{correct_principle}', which is described in option {correct_option_letter}. "
                f"The provided answer was {llm_answer}.")

# Execute the check and print the result
result = check_chemist_slang_answer()
print(result)