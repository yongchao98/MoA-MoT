def check_answer(llm_answer):
    """
    Checks the correctness of the answer to the analytical chemistry question.

    The function evaluates the given answer against the constraints provided in the question:
    1.  NMR: Two distinct peaks for the same proton.
    2.  LC: Two distinct peaks on a standard column.
    3.  MS: Both peaks show the same mass, consistent with the expected product.
    4.  Temperature: The two LC peaks are stable and resolved at elevated temperatures.
    """

    # Define the properties of each possible explanation
    options = {
        "A": {
            "name": "Enantiomers",
            "distinct_nmr_achiral": False,
            "distinct_lc_achiral": False,
            "same_mass": True,
            "stable_isomers": True
        },
        "B": {
            "name": "'Double coupling' product",
            "distinct_nmr_achiral": True,
            "distinct_lc_achiral": True,
            "same_mass": False, # This is a different molecule, so mass is different
            "stable_isomers": True
        },
        "C": {
            "name": "Precursor contaminant",
            "distinct_nmr_achiral": True,
            "distinct_lc_achiral": True,
            "same_mass": False, # This is a different molecule, so mass is different
            "stable_isomers": True
        },
        "D": {
            "name": "Diastereomers",
            "distinct_nmr_achiral": True,
            "distinct_lc_achiral": True,
            "same_mass": True,
            "stable_isomers": True # Not rotamers, which would coalesce at high temp
        }
    }

    # Extract the letter from the answer, e.g., "<<<D>>>" -> "D"
    try:
        answer_key = llm_answer.strip().replace("<", "").replace(">", "")
        if answer_key not in options:
            return f"Invalid answer format or option. The answer should be one of {list(options.keys())}."
    except Exception:
        return "Invalid answer format. Could not extract the answer key."

    chosen_option = options[answer_key]

    # Constraint 1: Two NMR peaks observed.
    # The explanation must account for two distinct signals in a standard (achiral) NMR.
    if not chosen_option["distinct_nmr_achiral"]:
        return (f"Incorrect. The answer '{chosen_option['name']}' is wrong because this species would not produce two distinct NMR peaks "
                "for the same proton in a standard achiral environment, but the question states two peaks were observed.")

    # Constraint 2: Two LC peaks observed.
    # The explanation must account for two separable peaks on a standard (achiral) LC column.
    if not chosen_option["distinct_lc_achiral"]:
        return (f"Incorrect. The answer '{chosen_option['name']}' is wrong because this species would not separate into two peaks "
                "on a standard achiral LC column, but the question states two peaks were observed.")

    # Constraint 3: Same mass spectrum for both peaks.
    # The explanation must involve species that have the same mass as the expected product.
    if not chosen_option["same_mass"]:
        return (f"Incorrect. The answer '{chosen_option['name']}' is wrong because this would be a different molecule with a different mass, "
                "which contradicts the LC-MS data showing both peaks have the same, correct mass.")

    # Constraint 4: Stable at elevated temperature.
    # This rules out rapidly interconverting species like rotamers, which would coalesce.
    # All given options are configurationally stable isomers or distinct molecules, so this check primarily reinforces the choice
    # of diastereomers over rotamers (which wasn't an option, but is a key chemical concept).
    if not chosen_option["stable_isomers"]:
         return (f"Incorrect. The answer '{chosen_option['name']}' is wrong because this phenomenon might not be stable at elevated temperatures, "
                 "contradicting the observation of two clearly defined LC peaks.")


    # If all checks pass, the answer is correct.
    # We can also verify that no other option works.
    correct_options = []
    for key, option in options.items():
        if (option["distinct_nmr_achiral"] and
            option["distinct_lc_achiral"] and
            option["same_mass"] and
            option["stable_isomers"]):
            correct_options.append(key)

    if answer_key in correct_options and len(correct_options) == 1:
        return "Correct"
    else:
        # This case should not be reached if the logic above is correct for a valid answer.
        return f"The answer {answer_key} is incorrect based on the analysis."


# The answer provided by the LLM
llm_answer_to_check = "<<<D>>>"

# Run the check
result = check_answer(llm_answer_to_check)
print(result)