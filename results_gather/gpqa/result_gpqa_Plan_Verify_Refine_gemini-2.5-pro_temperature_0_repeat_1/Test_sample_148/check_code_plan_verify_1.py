def check_chemistry_answer():
    """
    Checks the correctness of the answer to the peptide analysis question.

    This function simulates the reasoning process by encoding the properties
    of different chemical species and comparing them against the experimental
    observations provided in the question.
    """

    # --- 1. Define Experimental Observations from the Question ---
    observations = {
        "are_isomers": True,  # From LC-MS: two peaks with the same mass.
        "distinguishable_in_achiral_nmr": True,  # From NMR: two peaks for one alpha-proton.
        "separable_by_achiral_lc": True,  # From LC-MS: two clearly defined peaks.
    }

    # --- 2. Define Properties of Each Possible Explanation ---
    # These are based on fundamental principles of stereochemistry and analytical chemistry.
    properties = {
        'A': {  # 'Double coupling' product
            "are_isomers": False,  # Not an isomer, has a different mass.
            "distinguishable_in_achiral_nmr": True,
            "separable_by_achiral_lc": True,
            "reason_if_wrong": "A 'double coupled' product is not an isomer of the target molecule; it would have a different (larger) mass. This contradicts the MS data."
        },
        'B': {  # Precursor contamination
            "are_isomers": False,  # Not an isomer, has a different mass.
            "distinguishable_in_achiral_nmr": True,
            "separable_by_achiral_lc": True,
            "reason_if_wrong": "A precursor is not an isomer of the target molecule; it would have a different mass. This contradicts the MS data."
        },
        'C': {  # Diastereomers
            "are_isomers": True,   # Are isomers (same mass).
            "distinguishable_in_achiral_nmr": True,  # Have different chemical environments.
            "separable_by_achiral_lc": True,   # Have different physical properties.
            "reason_if_wrong": ""  # This is the correct answer.
        },
        'D': {  # Enantiomers
            "are_isomers": True,   # Are isomers (same mass).
            "distinguishable_in_achiral_nmr": False, # Indistinguishable in achiral solvent.
            "separable_by_achiral_lc": False,  # Not separable by standard (achiral) LC.
            "reason_if_wrong": "Enantiomers are not distinguishable by standard (achiral) NMR or separable by standard LC. This contradicts the observation of two distinct peaks in both analyses."
        }
    }

    # --- 3. Check the Provided Answer ---
    # The answer given by the LLM was 'C'.
    llm_answer = 'C'
    
    selected_option_properties = properties[llm_answer]
    
    # Compare the properties of the given answer with the observations.
    if (selected_option_properties["are_isomers"] == observations["are_isomers"] and
        selected_option_properties["distinguishable_in_achiral_nmr"] == observations["distinguishable_in_achiral_nmr"] and
        selected_option_properties["separable_by_achiral_lc"] == observations["separable_by_achiral_lc"]):
        
        # The answer 'C' is fully consistent with all observations.
        return "Correct"
    else:
        # If the answer were incorrect, determine the reason.
        mismatches = []
        if selected_option_properties["are_isomers"] != observations["are_isomers"]:
            mismatches.append("MS data (isomers vs. non-isomers)")
        if selected_option_properties["distinguishable_in_achiral_nmr"] != observations["distinguishable_in_achiral_nmr"]:
            mismatches.append("NMR data (distinguishable vs. indistinguishable)")
        if selected_option_properties["separable_by_achiral_lc"] != observations["separable_by_achiral_lc"]:
            mismatches.append("LC data (separable vs. non-separable)")
            
        return (f"Incorrect. The answer '{llm_answer}' is inconsistent with the following observations: {', '.join(mismatches)}. "
                f"Reason: {selected_option_properties['reason_if_wrong']}")

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)