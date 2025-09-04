def check_answer():
    """
    Checks the correctness of the provided answer by verifying chemical constraints.
    """
    # The final answer provided by the LLM to be checked.
    llm_answer = "C"

    # Define the properties of each candidate based on their IUPAC names.
    candidates = {
        "A": {
            "name": "3a,5,5-trimethyl-1,2,3,3a,5,6,7,8-octahydrocyclopenta[1,4]cyclobuta[1,2]benzene",
            "methyl_count": 3,
            "skeleton_rearranged": False
        },
        "B": {
            "name": "3a,4a,5,5-tetramethyl-2,3,3a,4,4a,5-hexahydro-1H-cyclobuta[1,2:1,4]di[5]annulene",
            "methyl_count": 4,
            "skeleton_rearranged": True
        },
        "C": {
            "name": "3a,4,5a-trimethyl-1,2,3,3a,5a,6,7,8-octahydrocyclopenta[c]pentalene",
            "methyl_count": 3,
            "skeleton_rearranged": True
        },
        "D": {
            "name": "3a,5-dimethyldecahydrocyclopenta[1,4]cyclobuta[1,2]benzene",
            "methyl_count": 2,
            "skeleton_rearranged": False
        }
    }

    # --- Constraint 1: Check the number of methyl groups ---
    # The final product must have 3 methyl groups (2 original + 1 from Wittig/TsOH steps).
    expected_methyl_count = 3
    if candidates[llm_answer]["methyl_count"] != expected_methyl_count:
        return (f"Incorrect. The answer {llm_answer} is wrong because it has "
                f"{candidates[llm_answer]['methyl_count']} methyl groups, but the final product must have "
                f"{expected_methyl_count} methyl groups.")

    # --- Constraint 2: Check for skeletal rearrangement ---
    # The final product must have a rearranged skeleton due to the relief of ring strain.
    if not candidates[llm_answer]["skeleton_rearranged"]:
        return (f"Incorrect. The answer {llm_answer} is wrong because it retains the original, strained "
                f"'cyclobuta' skeleton. The major product of the acid-catalyzed rearrangement should have a "
                f"rearranged, more stable skeleton.")

    # --- Verify that the chosen answer is the *only* correct one ---
    plausible_candidates = []
    for key, properties in candidates.items():
        if properties["methyl_count"] == expected_methyl_count and properties["skeleton_rearranged"]:
            plausible_candidates.append(key)
    
    if len(plausible_candidates) == 1 and plausible_candidates[0] == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. While answer {llm_answer} satisfies the basic constraints, the analysis is incomplete. "
                f"The set of plausible candidates is {plausible_candidates}, which should contain only one correct answer.")

# Run the check
result = check_answer()
print(result)