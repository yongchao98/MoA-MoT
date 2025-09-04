import re

def check_hplc_peaks():
    """
    This function simulates the chemical reasoning to determine the number of HPLC peaks
    and checks if the provided answer is correct.
    """

    # --- Step 1: Analyze the products of Reaction I ---
    # Starting Material: (S)-5-methoxyhexan-3-one. This is a chiral molecule with one existing stereocenter (C5).
    # Reaction: Reduction of the ketone at C3 creates a new stereocenter at C3.
    # Stereochemical Outcome: The original stereocenter (C5) is unaffected. The new stereocenter (C3) can be R or S.
    # This results in two products: (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol.
    # Relationship: These two products are diastereomers.
    products_reaction_I = {
        "description": "Two diastereomers",
        "total_stereoisomers": 2,
        "separable_groups_normal_phase": 2  # Diastereomers are separable
    }

    # --- Step 2: Analyze the products of Reaction II ---
    # Starting Material: Pentane-2,4-dione. This is an achiral molecule.
    # Reaction: Reduction of both ketones (C2 and C4) creates two new stereocenters.
    # Stereochemical Outcome: Since the starting material is achiral, all possible stereoisomers are formed.
    # Products: (2R, 4R)-pentane-2,4-diol, (2S, 4S)-pentane-2,4-diol, and (2R, 4S)-pentane-2,4-diol (meso).
    # Relationship: The (2R, 4R) and (2S, 4S) forms are a pair of enantiomers. The (2R, 4S) form is a meso compound, which is a diastereomer of the enantiomeric pair.
    # This gives a total of 3 distinct stereoisomers.
    products_reaction_II = {
        "description": "One pair of enantiomers and one meso compound",
        "total_stereoisomers": 3,
        "separable_groups_normal_phase": 2  # The enantiomeric pair co-elutes (1 peak), and the meso compound is separate (1 peak).
    }

    # --- Step 3: Calculate peaks for Normal-Phase HPLC ---
    # Principle: Separates diastereomers and constitutional isomers, but NOT enantiomers.
    # We sum the number of separable groups from each reaction, as the products of Rxn I and Rxn II are constitutionally different and will not co-elute.
    calculated_normal_phase_peaks = products_reaction_I["separable_groups_normal_phase"] + products_reaction_II["separable_groups_normal_phase"]

    # --- Step 4: Calculate peaks for Chiral HPLC ---
    # Principle: Separates all unique stereoisomers, including enantiomers.
    # We sum the total number of unique stereoisomers from each reaction.
    calculated_chiral_peaks = products_reaction_I["total_stereoisomers"] + products_reaction_II["total_stereoisomers"]

    # --- Step 5: Check the provided answer ---
    # The final answer provided by the LLM is D.
    final_answer_letter = "D"

    # Define the options from the question
    options = {
        "A": {"chiral": 3, "normal": 3},
        "B": {"chiral": 4, "normal": 2},
        "C": {"chiral": 3, "normal": 2},
        "D": {"chiral": 5, "normal": 4}
    }

    expected_peaks = options.get(final_answer_letter)

    if not expected_peaks:
        return f"Invalid answer letter '{final_answer_letter}' provided."

    # Compare calculated results with the expected answer
    errors = []
    if calculated_chiral_peaks != expected_peaks["chiral"]:
        errors.append(f"Chiral HPLC peak count is incorrect. Expected {expected_peaks['chiral']}, but calculated {calculated_chiral_peaks}.")
    
    if calculated_normal_phase_peaks != expected_peaks["normal"]:
        errors.append(f"Normal-Phase HPLC peak count is incorrect. Expected {expected_peaks['normal']}, but calculated {calculated_normal_phase_peaks}.")

    if not errors:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        reason = "The answer is incorrect for the following reason(s):\n"
        reason += "\n".join(errors)
        reason += "\n\n**Detailed Analysis Summary:**\n"
        reason += f"- **Reaction I** produces {products_reaction_I['total_stereoisomers']} diastereomers. This contributes {products_reaction_I['separable_groups_normal_phase']} peaks to normal-phase HPLC and {products_reaction_I['total_stereoisomers']} peaks to chiral HPLC.\n"
        reason += f"- **Reaction II** produces {products_reaction_II['total_stereoisomers']} total stereoisomers (a pair of enantiomers and a meso compound). This contributes {products_reaction_II['separable_groups_normal_phase']} peaks to normal-phase HPLC and {products_reaction_II['total_stereoisomers']} peaks to chiral HPLC.\n"
        reason += f"- **Total Normal-Phase Peaks:** {products_reaction_I['separable_groups_normal_phase']} + {products_reaction_II['separable_groups_normal_phase']} = {calculated_normal_phase_peaks}.\n"
        reason += f"- **Total Chiral Peaks:** {products_reaction_I['total_stereoisomers']} + {products_reaction_II['total_stereoisomers']} = {calculated_chiral_peaks}."
        return reason

# Run the check
result = check_hplc_peaks()
print(result)