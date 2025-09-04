def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by logically evaluating it against the problem's constraints.
    """
    llm_answer = 'C'

    # 1. Define the key experimental observations from the problem description.
    observations = {
        "two_nmr_peaks": True,        # Two peaks for the same alpha-proton.
        "two_lc_peaks": True,         # Two clearly defined peaks in LC.
        "mass_is_expected": True      # Both peaks have the same, expected mass.
    }

    # 2. Define the expected analytical properties for each possible answer.
    # This encodes the chemical principles.
    option_properties = {
        'A': {
            "description": "'Double coupling' product",
            "expected_mass": False,  # A double-coupled product would have a higher mass.
            "distinguishable_by_nmr": True,
            "separable_by_lc": True
        },
        'B': {
            "description": "Contamination with a precursor",
            "expected_mass": False,  # A precursor would have a different (likely lower) mass.
            "distinguishable_by_nmr": True,
            "separable_by_lc": True
        },
        'C': {
            "description": "Mixture of diastereoisomers",
            "expected_mass": True,   # Diastereomers are isomers and have the same mass.
            "distinguishable_by_nmr": True,  # Diastereomers have different NMR spectra.
            "separable_by_lc": True   # Diastereomers have different physical properties and are separable by standard LC.
        },
        'D': {
            "description": "Mixture of enantiomers",
            "expected_mass": True,   # Enantiomers are isomers and have the same mass.
            "distinguishable_by_nmr": False, # Enantiomers are indistinguishable by standard NMR (in achiral media).
            "separable_by_lc": False  # Enantiomers are not separable by standard (achiral) LC.
        }
    }

    # 3. Check the provided answer against the observations.
    chosen_option = option_properties.get(llm_answer)

    if not chosen_option:
        return f"Invalid answer option provided: {llm_answer}. Please choose from A, B, C, or D."

    # Constraint 1: Mass Spectrometry Data
    if chosen_option["expected_mass"] != observations["mass_is_expected"]:
        return (f"Incorrect. The answer '{llm_answer}' ({chosen_option['description']}) is inconsistent with the mass spectrometry data. "
                f"The MS shows the expected mass for both peaks, but a {chosen_option['description']} would have a different mass.")

    # Constraint 2: NMR Data
    if chosen_option["distinguishable_by_nmr"] != observations["two_nmr_peaks"]:
        return (f"Incorrect. The answer '{llm_answer}' ({chosen_option['description']}) is inconsistent with the NMR data. "
                f"The NMR shows two distinct peaks for one proton, but {chosen_option['description']} are indistinguishable by standard NMR.")

    # Constraint 3: LC Data
    if chosen_option["separable_by_lc"] != observations["two_lc_peaks"]:
        return (f"Incorrect. The answer '{llm_answer}' ({chosen_option['description']}) is inconsistent with the liquid chromatography data. "
                f"The LC shows two separate peaks, but {chosen_option['description']} are not separable by standard (achiral) LC.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_answer_correctness()
print(result)