def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by modeling the problem's constraints.

    The function defines the experimental observations as logical constraints and evaluates
    each multiple-choice option against them.
    """
    llm_answer = "C"

    # 1. Define the experimental observations from the question as constraints.
    # - LC-MS: Two peaks, same mass, mass is consistent with the expected molecule.
    # - NMR: Two peaks for the same proton, not due to coupling.
    # - LC: Two clearly defined peaks.
    constraints = {
        "are_isomers": True,  # From MS: same mass for both peaks
        "have_expected_mass": True,  # From MS: mass is consistent with the expected molecule
        "are_separable_by_standard_lc": True,  # From LC: two clearly defined peaks
        "are_distinguishable_by_standard_nmr": True  # From NMR: two peaks for the same proton
    }

    # 2. Define the properties of each possible explanation (the options A, B, C, D).
    explanations = {
        "A": {
            "name": "Mixture of enantiomers",
            "properties": {
                "are_isomers": True,
                "have_expected_mass": True,
                "are_separable_by_standard_lc": False,  # Enantiomers co-elute on achiral columns
                "are_distinguishable_by_standard_nmr": False  # Enantiomers are indistinguishable in achiral solvents
            }
        },
        "B": {
            "name": "Contaminated with a precursor",
            "properties": {
                "are_isomers": False,  # Different molecules
                "have_expected_mass": False, # Precursor has a different (lower) mass
                "are_separable_by_standard_lc": True,
                "are_distinguishable_by_standard_nmr": True
            }
        },
        "C": {
            "name": "Mixture of diastereoisomers",
            "properties": {
                "are_isomers": True,
                "have_expected_mass": True,
                "are_separable_by_standard_lc": True,  # Diastereomers have different physical properties
                "are_distinguishable_by_standard_nmr": True  # Diastereomers are chemically distinct
            }
        },
        "D": {
            "name": "'Double coupling' product",
            "properties": {
                "are_isomers": False,  # Different molecules
                "have_expected_mass": False, # Double coupling product has a different (higher) mass
                "are_separable_by_standard_lc": True,
                "are_distinguishable_by_standard_nmr": True
            }
        }
    }

    # 3. Determine the logically correct answer by finding which explanation matches all constraints.
    correct_option = None
    for option, data in explanations.items():
        if data["properties"] == constraints:
            correct_option = option
            break

    # 4. Compare the LLM's answer with the logically derived correct answer.
    if llm_answer == correct_option:
        return "Correct"
    else:
        # If the answer is wrong, explain which constraint was violated.
        chosen_explanation = explanations.get(llm_answer)
        if not chosen_explanation:
            return f"Incorrect. The answer '{llm_answer}' is not a valid option. The correct answer is '{correct_option}'."

        violated_constraints = []
        for key, value in constraints.items():
            if chosen_explanation["properties"][key] != value:
                violated_constraints.append(key)
        
        reason = f"The explanation '{chosen_explanation['name']}' is incorrect because it fails to satisfy the following observation(s): "
        
        if "are_isomers" in violated_constraints or "have_expected_mass" in violated_constraints:
            reason += "The mass spectrometry data shows both species are isomers with the expected mass. "
        if "are_separable_by_standard_lc" in violated_constraints:
            reason += "The liquid chromatography data shows two separable species. "
        if "are_distinguishable_by_standard_nmr" in violated_constraints:
            reason += "The NMR data shows two distinguishable species. "
            
        return f"Incorrect. {reason.strip()} The correct answer is '{correct_option}' ({explanations[correct_option]['name']}), which is consistent with all experimental data."

# Execute the check and print the result.
result = check_answer_correctness()
print(result)