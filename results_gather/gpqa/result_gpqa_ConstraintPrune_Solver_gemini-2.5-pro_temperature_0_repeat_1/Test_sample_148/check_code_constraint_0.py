def check_answer():
    """
    Checks the correctness of the answer by modeling the properties of each option
    and testing them against the experimental constraints from the question.
    """

    # Define the properties of each possible explanation based on chemical principles.
    options = {
        "A": {
            "name": "Enantiomers",
            "are_isomers": True,  # Isomers have the same mass.
            "distinguishable_by_achiral_nmr": False,  # Enantiomers have identical NMR spectra in achiral solvents.
            "separable_by_achiral_lc": False  # Enantiomers co-elute on standard (achiral) LC columns.
        },
        "B": {
            "name": "Precursor Contamination",
            "are_isomers": False, # A precursor is a different molecule with a different mass.
            "distinguishable_by_achiral_nmr": True,
            "separable_by_achiral_lc": True
        },
        "C": {
            "name": "'Double coupling' Product",
            "are_isomers": False, # A double-coupled product is a different molecule with a different (higher) mass.
            "distinguishable_by_achiral_nmr": True,
            "separable_by_achiral_lc": True
        },
        "D": {
            "name": "Diastereoisomers",
            "are_isomers": True,  # Diastereomers are isomers and have the same mass.
            "distinguishable_by_achiral_nmr": True,   # Diastereomers are chemically distinct and give different NMR spectra.
            "separable_by_achiral_lc": True   # Diastereomers have different physical properties and are separable by achiral LC.
        }
    }

    # The answer provided by the LLM.
    llm_answer = "D"

    # Define the constraints from the experimental observations.
    # Constraint 1: LC-MS shows two peaks with the same mass spectrum.
    # This implies the two compounds are isomers.
    constraint_ms = lambda props: props["are_isomers"]
    reason_ms = "The LC-MS data shows both peaks have the same mass, meaning they must be isomers. This option describes compounds with different masses."

    # Constraint 2: 1H NMR shows two distinct peaks for the same alpha-proton.
    # This implies the two compounds are distinguishable in a standard (achiral) NMR experiment.
    constraint_nmr = lambda props: props["distinguishable_by_achiral_nmr"]
    reason_nmr = "The NMR spectrum shows two distinct peaks, meaning the compounds must be distinguishable in an achiral environment. Enantiomers are not."

    # Constraint 3: LC analysis shows two clearly defined peaks.
    # This implies the two compounds are separable by standard (achiral) liquid chromatography.
    constraint_lc = lambda props: props["separable_by_achiral_lc"]
    reason_lc = "The LC analysis shows two separate peaks, meaning the compounds must be separable by a standard (achiral) column. Enantiomers are not."

    # Check if the LLM's chosen answer satisfies all constraints.
    selected_option_props = options.get(llm_answer)

    if not selected_option_props:
        return f"Invalid answer option '{llm_answer}' provided."

    if not constraint_ms(selected_option_props):
        return f"The answer '{llm_answer}' ({selected_option_props['name']}) is incorrect. It fails the Mass Spectrometry constraint. Reason: {reason_ms}"
    
    if not constraint_nmr(selected_option_props):
        return f"The answer '{llm_answer}' ({selected_option_props['name']}) is incorrect. It fails the NMR Spectroscopy constraint. Reason: {reason_nmr}"

    if not constraint_lc(selected_option_props):
        return f"The answer '{llm_answer}' ({selected_option_props['name']}) is incorrect. It fails the Liquid Chromatography constraint. Reason: {reason_lc}"

    # Verify that it is the *only* correct answer to ensure the logic is sound.
    correct_options = []
    for key, props in options.items():
        if constraint_ms(props) and constraint_nmr(props) and constraint_lc(props):
            correct_options.append(key)
    
    if len(correct_options) == 1 and correct_options[0] == llm_answer:
        return "Correct"
    else:
        # This case would be hit if the LLM's answer was correct but not unique, or if another answer was correct instead.
        return f"The provided answer '{llm_answer}' is incorrect. The only option that satisfies all experimental constraints is '{correct_options[0]}'."

# Execute the check and print the result.
result = check_answer()
print(result)