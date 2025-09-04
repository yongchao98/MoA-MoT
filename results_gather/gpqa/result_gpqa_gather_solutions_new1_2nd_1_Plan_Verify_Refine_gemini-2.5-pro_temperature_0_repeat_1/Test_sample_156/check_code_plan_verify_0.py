def check_correctness_of_retrovirus_kit_answer():
    """
    This function checks the correctness of the selected answer for designing a
    molecular diagnostic kit for a retrovirus.

    It evaluates each option based on three key constraints from the question:
    1. Must be a 'molecular' test (detects nucleic acids).
    2. Must be biologically correct for a 'retrovirus' (RNA genome).
    3. Must be 'quick and accurate' for an outbreak scenario.
    """

    # The final answer provided by the LLM that we need to check.
    # Note: The provided text has some answers as A, B, C, and some as D.
    # The final consolidated answer at the end of the prompt is D.
    llm_answer = 'D'

    # Define the properties of each option based on scientific principles.
    options = {
        'A': {
            'is_molecular': False,
            'reason_molecular': "Fails constraint: It describes an ELISA kit, which is an immunological test detecting antibodies (proteins), not a molecular test.",
            'is_biologically_correct': False,
            'reason_biology': "Fails constraint: It detects the host's immune response (IgG), not the virus's genetic material directly.",
            'is_quick_and_accurate': False,
            'reason_quick': "Fails constraint: IgG antibodies are markers of a late-stage or past infection, making the test unsuitable for 'quick' early diagnosis."
        },
        'B': {
            'is_molecular': True,
            'is_biologically_correct': False,
            'reason_biology': "Fails constraint: It proposes direct 'DNA sequencing' on a retrovirus, which is biologically incorrect as retroviruses have an RNA genome. It misses the essential reverse transcription step.",
            'is_quick_and_accurate': True, # PCR itself is quick, but the premise is flawed.
        },
        'C': {
            'is_molecular': True,
            'is_biologically_correct': False,
            'reason_biology': "Fails constraint: The identification method is scientifically invalid. One cannot determine a virus's precise genetic sequence from non-specific symptoms to design a PCR test.",
            'is_quick_and_accurate': False, # Nested PCR is slower and more contamination-prone than real-time PCR.
            'reason_quick': "Fails constraint: Nested PCR is less suitable for a 'quick' high-throughput kit compared to real-time PCR."
        },
        'D': {
            'is_molecular': True,
            'is_biologically_correct': True,
            'is_quick_and_accurate': True,
        }
    }

    # Find the truly correct option by checking all constraints.
    scientifically_correct_option = None
    for option_id, properties in options.items():
        if all(properties.get(key, True) for key in ['is_molecular', 'is_biologically_correct', 'is_quick_and_accurate']):
            scientifically_correct_option = option_id
            break

    # Check if the LLM's answer matches the scientifically correct one.
    if llm_answer == scientifically_correct_option:
        return "Correct"
    else:
        # Provide a reason why the LLM's answer is wrong.
        if llm_answer not in options:
            return f"Invalid answer choice '{llm_answer}' provided."

        chosen_option_properties = options[llm_answer]
        reasons_for_failure = []
        if not chosen_option_properties.get('is_molecular', True):
            reasons_for_failure.append(chosen_option_properties['reason_molecular'])
        if not chosen_option_properties.get('is_biologically_correct', True):
            reasons_for_failure.append(chosen_option_properties['reason_biology'])
        if not chosen_option_properties.get('is_quick_and_accurate', True):
            reasons_for_failure.append(chosen_option_properties['reason_quick'])

        if not reasons_for_failure:
             return f"Incorrect. The provided answer {llm_answer} is wrong, but the checker could not find a specific reason. The correct answer is {scientifically_correct_option}."

        return f"Incorrect. The provided answer '{llm_answer}' is wrong because it violates key constraints of the question. Specifically: {' '.join(reasons_for_failure)}"

# Run the check and print the result.
print(check_correctness_of_retrovirus_kit_answer())