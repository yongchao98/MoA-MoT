def check_chemistry_hplc_answer():
    """
    This function programmatically verifies the answer to the organic chemistry HPLC problem.
    It models the products of each reaction and simulates the separation on both
    normal-phase and chiral HPLC columns.
    """

    # --- Step 1: Define the products from Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one (chiral)
    # Reaction: Reduction of ketone at C3 creates a new stereocenter.
    # Products: A pair of diastereomers, (3R, 5S)- and (3S, 5S)-5-methoxyhexan-3-ol.
    products_rxn1 = {
        "description": "A pair of diastereomers",
        "num_stereoisomers": 2,
        "num_normal_hplc_peaks": 2  # Diastereomers are separable by normal HPLC
    }

    # --- Step 2: Define the products from Reaction II ---
    # Reactant: Pentane-2,4-dione (achiral)
    # Reaction: Reduction of both ketones creates two new stereocenters.
    # Products: A pair of enantiomers ((2R,4R) and (2S,4S)) and a meso compound ((2R,4S)).
    products_rxn2 = {
        "description": "An enantiomeric pair and a meso compound",
        "num_stereoisomers": 3,
        "num_normal_hplc_peaks": 2  # Enantiomers co-elute (1 peak), meso is a diastereomer (1 peak)
    }

    # --- Step 3: Calculate the total number of peaks for each HPLC type ---
    # The products of Reaction I and Reaction II are constitutionally different
    # (different molecular formulas, C7H16O2 vs C5H12O2) and will not co-elute.
    # Therefore, we can sum the peaks from each reaction.

    # Chiral HPLC separates all unique stereoisomers.
    total_chiral_peaks = products_rxn1["num_stereoisomers"] + products_rxn2["num_stereoisomers"]

    # Normal-phase HPLC separates diastereomers and constitutional isomers, but not enantiomers.
    total_normal_peaks = products_rxn1["num_normal_hplc_peaks"] + products_rxn2["num_normal_hplc_peaks"]

    # --- Step 4: Define the options and the provided answer ---
    options = {
        "A": {"chiral": 5, "normal": 4},
        "B": {"chiral": 3, "normal": 2},
        "C": {"chiral": 3, "normal": 3},
        "D": {"chiral": 4, "normal": 2}
    }
    llm_answer_choice = "A"

    # --- Step 5: Check if the calculated results match the provided answer ---
    expected_chiral = options[llm_answer_choice]["chiral"]
    expected_normal = options[llm_answer_choice]["normal"]

    if total_chiral_peaks == expected_chiral and total_normal_peaks == expected_normal:
        return "Correct"
    else:
        error_messages = []
        if total_chiral_peaks != expected_chiral:
            error_messages.append(f"Chiral HPLC peak count is wrong. Calculated {total_chiral_peaks}, but answer {llm_answer_choice} implies {expected_chiral}.")
        if total_normal_peaks != expected_normal:
            error_messages.append(f"Normal-phase HPLC peak count is wrong. Calculated {total_normal_peaks}, but answer {llm_answer_choice} implies {expected_normal}.")
        return "Incorrect. " + " ".join(error_messages)

# Execute the check
result = check_chemistry_hplc_answer()
print(result)