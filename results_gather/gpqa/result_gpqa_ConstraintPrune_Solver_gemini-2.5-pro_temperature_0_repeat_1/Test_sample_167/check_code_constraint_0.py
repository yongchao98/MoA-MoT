def check_genomics_error_answer():
    """
    Checks the correctness of the answer by evaluating which genomics errors
    are considered "difficult-to-spot" (i.e., cause silent failures).
    """
    # The given answer from the LLM
    llm_answer = "A"

    # Define the issues and their classification based on common bioinformatics knowledge.
    # True = difficult to spot (silent failure), False = easy to spot (loud failure).
    issue_classification = {
        1: {"name": "Mutually incompatible data formats", "is_difficult": False},
        2: {"name": "The 'chr' / 'no chr' confusion", "is_difficult": True},
        3: {"name": "Reference assembly mismatch", "is_difficult": True},
        4: {"name": "Incorrect ID conversion", "is_difficult": True}
    }

    # Determine the correct set of issues based on the "difficult-to-spot" criterion.
    correct_issue_set = {
        issue_id for issue_id, props in issue_classification.items() if props["is_difficult"]
    }
    # Expected correct_issue_set: {2, 3, 4}

    # Define the sets corresponding to the multiple-choice options.
    options = {
        "A": {2, 3, 4},
        "B": {3, 4},
        "C": {2, 3},
        "D": {1, 2, 3, 4}
    }

    # Check if the provided answer key exists in the options.
    if llm_answer not in options:
        return f"Incorrect. The answer '{llm_answer}' is not a valid option."

    # Get the set of issues for the given answer.
    llm_answer_set = options[llm_answer]

    # Compare the answer's set with the correct set.
    if llm_answer_set == correct_issue_set:
        return "Correct"
    else:
        # Identify which issues were incorrectly included or excluded.
        incorrectly_included = llm_answer_set - correct_issue_set
        incorrectly_excluded = correct_issue_set - llm_answer_set

        reasons = []
        if incorrectly_included:
            for issue_id in incorrectly_included:
                reasons.append(
                    f"Issue {issue_id} ('{issue_classification[issue_id]['name']}') was incorrectly included. "
                    "This type of error typically causes an immediate and obvious failure (e.g., a program crash), "
                    "making it easy to spot, not difficult."
                )
        
        if incorrectly_excluded:
            for issue_id in incorrectly_excluded:
                 reasons.append(
                    f"Issue {issue_id} ('{issue_classification[issue_id]['name']}') was incorrectly excluded. "
                    "This is a classic source of difficult-to-spot errors, as the analysis may run without "
                    "crashing but produce scientifically invalid results."
                )

        return f"Incorrect. The answer '{llm_answer}' is wrong for the following reason(s): {' '.join(reasons)}"

# Run the check
result = check_genomics_error_answer()
print(result)