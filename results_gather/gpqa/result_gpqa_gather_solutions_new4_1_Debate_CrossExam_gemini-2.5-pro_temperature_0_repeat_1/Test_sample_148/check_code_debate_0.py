def check_answer_correctness():
    """
    Checks the correctness of the final answer for the peptide analysis problem.

    The function codifies the experimental observations as logical constraints
    and evaluates each possible explanation against these constraints.
    """

    # --- Step 1: Define the final answer to be checked ---
    # The provided final answer is <<<B>>>, which corresponds to a mixture of diastereoisomers.
    final_answer = "B"

    # --- Step 2: Define the experimental observations from the question ---
    # These are the key facts from the NMR and LC-MS data.
    observations = {
        "nmr_shows_multiple_signals_for_one_proton": True,
        "lc_separates_into_multiple_peaks": True,
        "ms_shows_same_mass_as_expected_product": True
    }

    # --- Step 3: Define the properties of each possible explanation based on chemical principles ---
    # This maps each option to its expected analytical signature.
    explanations = {
        "A": {  # Contamination with a precursor
            "description": "Contamination with a precursor",
            "has_same_mass_as_product": False,
            "is_distinguishable_by_achiral_nmr_lc": True
        },
        "B": {  # Mixture of diastereoisomers
            "description": "Mixture of diastereoisomers",
            "has_same_mass_as_product": True,
            "is_distinguishable_by_achiral_nmr_lc": True
        },
        "C": {  # 'Double coupling' side reaction
            "description": "'Double coupling' side reaction",
            "has_same_mass_as_product": False,
            "is_distinguishable_by_achiral_nmr_lc": True
        },
        "D": {  # Mixture of enantiomers
            "description": "Mixture of enantiomers",
            "has_same_mass_as_product": True,
            "is_distinguishable_by_achiral_nmr_lc": False
        }
    }

    # --- Step 4: Determine the logically correct option based on observations ---
    correct_option = None
    for option, properties in explanations.items():
        # Check if the explanation matches all observations
        if (properties["has_same_mass_as_product"] == observations["ms_shows_same_mass_as_expected_product"] and
            properties["is_distinguishable_by_achiral_nmr_lc"] == observations["nmr_shows_multiple_signals_for_one_proton"]):
            correct_option = option
            break

    # --- Step 5: Compare the final answer with the logically derived correct option ---
    if final_answer == correct_option:
        return "Correct"
    else:
        # If the answer is incorrect, provide a detailed reason.
        if final_answer not in explanations:
            return f"The provided answer '{final_answer}' is not one of the valid options (A, B, C, D)."

        chosen_explanation = explanations[final_answer]
        
        # Check which constraint the incorrect answer violates.
        if chosen_explanation["has_same_mass_as_product"] != observations["ms_shows_same_mass_as_expected_product"]:
            return (f"The answer '{final_answer}' ({chosen_explanation['description']}) is incorrect. "
                    f"This explanation is inconsistent with the mass spectrometry data. The data shows both peaks have the same mass as the expected product, "
                    f"whereas a {chosen_explanation['description'].lower()} would have a different mass.")

        if chosen_explanation["is_distinguishable_by_achiral_nmr_lc"] != observations["nmr_shows_multiple_signals_for_one_proton"]:
            return (f"The answer '{final_answer}' ({chosen_explanation['description']}) is incorrect. "
                    f"This explanation is inconsistent with the NMR and LC data. The data shows two distinct, separable species, "
                    f"but {chosen_explanation['description'].lower()} are not distinguishable by standard (achiral) NMR and LC techniques.")
        
        return f"The answer '{final_answer}' is incorrect for the reasons outlined."

# Execute the check and print the result.
result = check_answer_correctness()
print(result)