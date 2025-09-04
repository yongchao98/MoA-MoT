def check_answer_correctness():
    """
    Checks the correctness of the selected answer for designing a retrovirus diagnostic kit.
    """
    question_constraints = {
        "kit_type": "molecular",  # Must detect nucleic acids
        "pathogen_type": "retrovirus",  # RNA genome, requires reverse transcription
        "detection_speed": "quick"  # Must be suitable for early and rapid diagnosis
    }

    # Properties of each option based on the question's text
    options_analysis = {
        "A": {
            "is_molecular": True,
            "handles_retrovirus_correctly": False,
            "is_quick": False,
            "reason": "The initial identification step is scientifically invalid; you cannot design a PCR kit from symptoms. Nested PCR is also generally slower and more prone to contamination than real-time PCR."
        },
        "B": {
            "is_molecular": True,
            "handles_retrovirus_correctly": False,
            "is_quick": True,
            "reason": "It incorrectly assumes the retrovirus has a DNA genome. It misses the essential reverse transcription step (RNA to cDNA) required before sequencing or PCR."
        },
        "C": {
            "is_molecular": True,
            "handles_retrovirus_correctly": True,
            "is_quick": True,
            "reason": "This option correctly describes the standard workflow: convert RNA to cDNA for sequencing, then use the sequence to design a rapid and sensitive real-time PCR test."
        },
        "D": {
            "is_molecular": False,
            "handles_retrovirus_correctly": False, # Does not target nucleic acid
            "is_quick": False,
            "reason": "This describes a serological (immunological) test, not a molecular one. Furthermore, detecting IgG antibodies is not suitable for quick, early diagnosis as they appear late in the infection."
        }
    }

    # The answer provided by the LLM being checked
    llm_answer = "C"

    selected_option = options_analysis.get(llm_answer)

    if not selected_option:
        return f"Invalid answer option: {llm_answer}"

    # Check against constraints
    if not selected_option["is_molecular"]:
        return f"Incorrect. The question requires a '{question_constraints['kit_type']}' kit. Option {llm_answer}'s method is not molecular. Reason: {selected_option['reason']}"

    if not selected_option["handles_retrovirus_correctly"]:
        return f"Incorrect. The method in option {llm_answer} does not correctly handle a '{question_constraints['pathogen_type']}'. Reason: {selected_option['reason']}"

    if not selected_option["is_quick"]:
        return f"Incorrect. The question requires '{question_constraints['detection_speed']}' detection. Option {llm_answer} is not suitable for this. Reason: {selected_option['reason']}"

    # If all constraints are met
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)