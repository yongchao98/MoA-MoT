def check_answer():
    """
    This function checks the correctness of the given answer by logically evaluating
    the observations from the problem against the properties of each possible explanation.
    """
    llm_answer = 'D'

    # --- Observations from the problem statement ---
    # 1. NMR: Two peaks for the same alpha-proton, not due to spin-spin coupling.
    # This implies two distinct chemical species with different magnetic environments.
    has_distinguishable_nmr_spectra = True

    # 2. LC-MS: Two clearly defined peaks.
    # This implies two physically separable species.
    is_separable_by_lc = True

    # 3. LC-MS: Both peaks have the same mass spectrum.
    # This implies the two species are isomers (same mass).
    are_isomers = True

    # --- Evaluate each option based on chemical principles ---

    # A) Mixture of enantiomers
    # Enantiomers are NOT separable by standard LC and are NOT distinguishable by standard NMR.
    # This contradicts the observations.
    a_is_correct = (not is_separable_by_lc) and (not has_distinguishable_nmr_spectra) and are_isomers

    # B) Contaminated with a precursor
    # A precursor is a different molecule, not an isomer. It would have a different mass.
    # This contradicts the "same mass" observation.
    b_is_correct = not are_isomers

    # C) 'Double coupling' has occurred
    # This would create a side-product with a different chemical formula and mass, not an isomer.
    # This contradicts the "same mass" observation.
    c_is_correct = not are_isomers

    # D) Mixture of diastereoisomers
    # Diastereomers are isomers (same mass), are separable by LC, and have distinct NMR spectra.
    # This matches all observations.
    d_is_correct = is_separable_by_lc and has_distinguishable_nmr_spectra and are_isomers

    # --- Determine the logically correct answer ---
    correct_option = None
    if d_is_correct:
        correct_option = 'D'
    elif a_is_correct:
        correct_option = 'A'
    elif b_is_correct:
        correct_option = 'B'
    elif c_is_correct:
        correct_option = 'C'

    # --- Compare the LLM's answer with the logically derived answer ---
    if llm_answer == correct_option:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is '{llm_answer}', but the correct answer is '{correct_option}'.\n"
                f"Reasoning:\n"
                f"The observations are: separable by LC ({is_separable_by_lc}), distinguishable by NMR ({has_distinguishable_nmr_spectra}), and are isomers ({are_isomers}).\n"
                f"- Option A (Enantiomers) fails because they are not separable by standard LC or distinguishable by NMR.\n"
                f"- Options B (Precursor) and C (Side-product) fail because they are not isomers and would have a different mass.\n"
                f"- Option D (Diastereomers) is the only choice that fits all three observations.")

# Execute the check
result = check_answer()
print(result)