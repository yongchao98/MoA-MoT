def check_correctness():
    """
    Checks the correctness of the LLM's answer by evaluating it against
    the constraints of the chemistry problem.
    """
    llm_answer = "B"
    options = {
        "A": "The compounds they are working with have similar optical rotations.",
        "B": "The compounds they are working with have similar polarities.",
        "C": "The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.",
        "D": "The compounds they are working with have similar boiling points."
    }

    # --- Knowledge Base ---
    # We define the properties of each option based on chemical principles.
    # is_separation_principle: Is this property the basis for a common lab separation?
    # vernacular_match: Does the phrase "on top of each other" commonly describe a failure of the associated technique?
    analysis = {
        "A": {
            "is_separation_principle": False,
            "reason_separation": "Optical rotation is a characterization property, not a principle for separating bulk compounds.",
            "vernacular_match": False
        },
        "B": {
            "is_separation_principle": True,
            "technique": "Chromatography",
            "vernacular_match": True,
            "reason_vernacular": "In chromatography (TLC, column), compounds with similar polarities co-elute, resulting in spots/bands that are visually 'on top of each other'."
        },
        "C": {
            "is_separation_principle": False,
            "reason_separation": "While interactions are part of separation, 'bonding to each other' is not the primary principle used for separation. The key is the differential interaction with a stationary/mobile phase.",
            "vernacular_match": False
        },
        "D": {
            "is_separation_principle": True,
            "technique": "Distillation",
            "vernacular_match": False,
            "reason_vernacular": "A failed distillation (due to similar boiling points) is described as 'co-distilling' or having mixed fractions, not as compounds being 'on top of each other'."
        }
    }

    # --- Evaluation ---
    correct_option = None
    for key, props in analysis.items():
        if props["is_separation_principle"] and props["vernacular_match"]:
            correct_option = key
            break

    if correct_option is None:
        # This case should not be reached if the logic is sound for this question.
        return "Error in checking logic: No option satisfies all constraints."

    if llm_answer == correct_option:
        return "Correct"
    else:
        failure_reason = f"The provided answer '{llm_answer}' is incorrect. The correct answer is '{correct_option}'.\n\n"
        failure_reason += f"Here's the breakdown:\n"
        failure_reason += f"1. The problem describes a failed separation. This eliminates option A (characterization) and C (not a primary separation principle).\n"
        failure_reason += f"2. We are left with B (polarity/chromatography) and D (boiling point/distillation).\n"
        failure_reason += f"3. The key phrase 'on top of each other' is specific jargon for chromatography, where unresolved compounds create overlapping visual spots or bands.\n"
        failure_reason += f"4. This phrase is not used for distillation. Therefore, the compounds must have similar polarities (B), leading to a failed chromatographic separation.\n"
        failure_reason += f"The answer '{llm_answer}' does not fit this reasoning."
        return failure_reason

# Execute the check and print the result.
result = check_correctness()
print(result)