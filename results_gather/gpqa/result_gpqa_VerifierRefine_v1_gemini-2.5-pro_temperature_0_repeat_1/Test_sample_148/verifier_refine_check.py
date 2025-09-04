def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer by logically evaluating each option
    against the experimental observations provided in the question.
    """

    # 1. Define the experimental observations from the question
    observations = {
        "nmr_peaks_distinguishable": True,  # Two peaks for the same alpha-proton
        "lc_peaks_separable": True,         # Two clearly defined LC peaks
        "mass_of_species": "same_as_expected", # Both peaks have the same mass, consistent with the expected molecule
        "analysis_conditions": "achiral"    # Implied by standard NMR and LC-MS unless specified otherwise
    }

    # 2. Define the properties of the chemical species described in each option
    options_properties = {
        "A": {
            "description": "Contaminated with a precursor",
            "mass": "different_from_expected",
            "distinguishable_in_achiral_nmr": True,
            "separable_by_achiral_lc": True,
        },
        "B": {
            "description": "Mixture of enantiomers",
            "mass": "same_as_expected",
            "distinguishable_in_achiral_nmr": False,
            "separable_by_achiral_lc": False,
        },
        "C": {
            "description": "'Double coupling' byproduct",
            "mass": "different_from_expected", # Specifically, a higher mass
            "distinguishable_in_achiral_nmr": True,
            "separable_by_achiral_lc": True,
        },
        "D": {
            "description": "Mixture of diastereoisomers",
            "mass": "same_as_expected",
            "distinguishable_in_achiral_nmr": True,
            "separable_by_achiral_lc": True,
        }
    }

    # The answer provided by the LLM to be checked
    llm_answer = "D"

    # 3. Perform the logical validation
    consistent_options = []
    reasons_for_elimination = {}

    for option, properties in options_properties.items():
        is_consistent = True
        inconsistencies = []

        # Check 1: Mass consistency
        if properties["mass"] != observations["mass_of_species"]:
            is_consistent = False
            inconsistencies.append("it contradicts the mass spectrometry data, which shows both species have the same mass as the expected molecule.")

        # Check 2: NMR distinguishability under achiral conditions
        if properties["distinguishable_in_achiral_nmr"] != observations["nmr_peaks_distinguishable"]:
            is_consistent = False
            inconsistencies.append("it contradicts the NMR data. Enantiomers (Option B) would not produce two distinct peaks in a standard (achiral) NMR.")

        # Check 3: LC separability under achiral conditions
        if properties["separable_by_achiral_lc"] != observations["lc_peaks_separable"]:
            is_consistent = False
            inconsistencies.append("it contradicts the LC data. Enantiomers (Option B) are not separable on a standard (achiral) LC column.")

        if is_consistent:
            consistent_options.append(option)
        else:
            reasons_for_elimination[option] = " ".join(inconsistencies)

    # 4. Formulate the final verdict
    if llm_answer in consistent_options:
        if len(consistent_options) == 1:
            # The LLM's answer is the only one that fits all criteria.
            return "Correct"
        else:
            # This case indicates an ambiguous question, but not for this problem.
            return f"Incorrect. The answer {llm_answer} is plausible, but other options {consistent_options} are also consistent with the data."
    else:
        # The LLM's answer is incorrect. Provide the reason.
        reason = reasons_for_elimination.get(llm_answer, "The provided answer is not a valid option.")
        return f"Incorrect. The answer {llm_answer} is wrong because {reason}"

# Execute the check and print the result
result = check_answer_correctness()
print(result)