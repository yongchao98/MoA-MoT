def check_answer():
    """
    This function checks the correctness of the final answer by simulating the logical deduction process
    based on the experimental evidence provided in the question.
    """

    # 1. Define the properties of each possible explanation based on chemical principles.
    # 'is_isomer': Does it have the same mass as the target molecule?
    # 'is_distinguishable': Is it distinguishable by standard (achiral) NMR and separable by standard LC?
    properties = {
        'precursor': {
            'is_isomer': False,  # Has a different (lower) mass
            'is_distinguishable': True
        },
        'enantiomers': {
            'is_isomer': True,   # Has the same mass
            'is_distinguishable': False # Not distinguishable in an achiral environment
        },
        'diastereoisomers': {
            'is_isomer': True,   # Has the same mass
            'is_distinguishable': True
        },
        'double_coupling': {
            'is_isomer': False, # Has a different (higher) mass
            'is_distinguishable': True
        }
    }

    # 2. Define the experimental observations from the question as constraints.
    # Constraint from MS: "Both peaks have the same mass spectrum" -> The species are isomers.
    constraint_ms = {'is_isomer': True}
    # Constraint from NMR/LC: "two peaks" in both -> The species are distinguishable/separable.
    constraint_nmr_lc = {'is_distinguishable': True}

    # 3. Map the question's options to our defined explanations.
    # A) precursor, B) enantiomers, C) diastereoisomers, D) double coupling
    question_options = {
        'A': 'precursor',
        'B': 'enantiomers',
        'C': 'diastereoisomers',
        'D': 'double_coupling'
    }

    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'C'

    # 4. Perform the logical deduction.
    
    # Start with all possible options.
    possible_options = set(question_options.keys())
    
    # Apply Constraint 1 (from Mass Spectrometry)
    options_after_ms = set()
    for option in possible_options:
        explanation_type = question_options[option]
        if properties[explanation_type]['is_isomer'] == constraint_ms['is_isomer']:
            options_after_ms.add(option)
    
    # Apply Constraint 2 (from NMR and LC)
    options_after_nmr_lc = set()
    for option in options_after_ms:
        explanation_type = question_options[option]
        if properties[explanation_type]['is_distinguishable'] == constraint_nmr_lc['is_distinguishable']:
            options_after_nmr_lc.add(option)

    # The set should now contain only the single correct answer.
    if len(options_after_nmr_lc) != 1:
        return f"Logic Error: The analysis resulted in {len(options_after_nmr_lc)} possible answers, but there should be exactly one."

    correct_option = options_after_nmr_lc.pop()

    # 5. Compare the deduced correct answer with the LLM's final answer.
    if llm_final_answer == correct_option:
        return "Correct"
    else:
        # Provide a reason why the LLM's answer is wrong.
        llm_explanation = question_options.get(llm_final_answer, "an invalid option")
        correct_explanation = question_options[correct_option]
        
        reason = f"Incorrect. The provided answer is '{llm_final_answer}' ({llm_explanation}), but the correct answer is '{correct_option}' ({correct_explanation}).\n"
        
        # Check which constraint the LLM's answer failed.
        if not properties[llm_explanation]['is_isomer']:
            reason += f"Reason: A {llm_explanation} is not an isomer (it has a different mass), but the MS data showed both compounds have the same mass."
        elif not properties[llm_explanation]['is_distinguishable']:
            reason += f"Reason: {llm_explanation.capitalize()} are not distinguishable by standard NMR or separable by standard LC, but the experimental data showed two distinct peaks in both analyses."
        
        return reason

# Execute the check
result = check_answer()
print(result)