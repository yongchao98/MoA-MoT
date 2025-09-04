def check_genomics_answer():
    """
    Checks the correctness of the answer to the genomics data analysis question
    by encoding expert knowledge as a ground truth.
    """
    # The question asks for issues that are BOTH "common" AND "difficult-to-spot".
    # "Difficult-to-spot" means the error can lead to silent failures (plausible but wrong results)
    # rather than just obvious crashes.

    # Define the issues and their properties based on established knowledge in bioinformatics.
    # This represents the ground truth for the check.
    issues_ground_truth = {
        1: {
            "name": "Mutually incompatible data formats",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reason": "Includes subtle but common errors like 0-based vs. 1-based coordinate confusion, which cause silent off-by-one errors."
        },
        2: {
            "name": "The 'chr' / 'no chr' confusion",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reason": "A classic silent failure where tools find no overlaps and return a 'zero' result, which can be misinterpreted as a valid biological finding."
        },
        3: {
            "name": "Reference assembly mismatch",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reason": "An insidious error where analysis completes without technical warnings, but results are biologically invalid due to incorrect genomic coordinates."
        },
        4: {
            "name": "Incorrect ID conversion",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reason": "Error happens silently (e.g., in spreadsheets) and corrupts input data, leading to flawed downstream analysis on a correct-looking dataset."
        }
    }

    # Determine the set of correct issue indices based on the ground truth criteria.
    correct_issue_indices = {
        idx for idx, properties in issues_ground_truth.items()
        if properties["is_common"] and properties["is_difficult_to_spot"]
    }

    # The options as presented in the final answer's context.
    options_map = {
        "A": {3, 4},
        "B": {1, 2, 3, 4},  # "All of the above"
        "C": {2, 3},
        "D": {2, 3, 4}
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_choice = "B"

    # Get the set of issues corresponding to the LLM's answer.
    llm_answer_issues = options_map.get(llm_answer_choice)

    # Compare the LLM's answer with the ground truth.
    if llm_answer_issues == correct_issue_indices:
        return "Correct"
    else:
        missing_issues = correct_issue_indices - llm_answer_issues
        extra_issues = llm_answer_issues - correct_issue_indices

        reason = f"The answer '{llm_answer_choice}' is incorrect. "
        
        if missing_issues:
            missing_names = [issues_ground_truth[i]['name'] for i in sorted(list(missing_issues))]
            reason += f"It incorrectly omits the following issues which are common and difficult-to-spot: {', '.join(missing_names)}. "
        
        if extra_issues:
            extra_names = [issues_ground_truth[i]['name'] for i in sorted(list(extra_issues))]
            reason += f"It incorrectly includes issues that are not considered both common and difficult-to-spot: {', '.join(extra_names)}. "
            
        reason += "The correct answer should include all four issues because each is a well-documented, common source of errors that can 'fail silently' without causing an obvious program crash."
        return reason.strip()

# Execute the check
result = check_genomics_answer()
print(result)