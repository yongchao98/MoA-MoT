def check_genomics_question_answer():
    """
    Checks the correctness of the answer to the genomics data analysis question.

    The function codifies the expert knowledge presented in the provided answer's reasoning:
    1. It defines the properties of each potential issue (is it common? is it difficult-to-spot?).
    2. It identifies the set of issues that meet the question's criteria.
    3. It compares this "correct" set with the set represented by the given answer choice.
    """

    # Step 1: Define the issues and the expert consensus based on the provided reasoning.
    # The reasoning concludes that all four issues are both common and difficult-to-spot.
    issues_properties = {
        1: {"name": "Mutually incompatible data formats", "is_common": True, "is_difficult_to_spot": True},
        2: {"name": "The 'chr' / 'no chr' confusion", "is_common": True, "is_difficult_to_spot": True},
        3: {"name": "Reference assembly mismatch", "is_common": True, "is_difficult_to_spot": True},
        4: {"name": "Incorrect ID conversion", "is_common": True, "is_difficult_to_spot": True}
    }

    # Step 2: Define the mapping from the multiple-choice options to the issue numbers.
    options_map = {
        "A": {1, 2, 3, 4},  # All of the above
        "B": {3, 4},
        "C": {2, 3},
        "D": {2, 3, 4}
    }

    # The final answer provided by the LLM to be checked.
    given_answer_key = "A"

    # Step 3: Determine the correct set of issues based on the defined properties.
    # The question asks for issues that are BOTH common AND difficult-to-spot.
    correct_issue_numbers = set()
    for issue_num, properties in issues_properties.items():
        if properties["is_common"] and properties["is_difficult_to_spot"]:
            correct_issue_numbers.add(issue_num)

    # Step 4: Get the set of issues corresponding to the given answer.
    if given_answer_key not in options_map:
        return f"Incorrect. The provided answer key '{given_answer_key}' is not a valid option."

    answered_issue_numbers = options_map[given_answer_key]

    # Step 5: Compare the correct set with the answered set.
    if correct_issue_numbers == answered_issue_numbers:
        return "Correct"
    else:
        # Provide a detailed reason for the mismatch.
        missing_issues = correct_issue_numbers - answered_issue_numbers
        extra_issues = answered_issue_numbers - correct_issue_numbers

        error_messages = []
        if missing_issues:
            missing_names = [f"'{issues_properties[i]['name']}' (Issue {i})" for i in sorted(list(missing_issues))]
            error_messages.append(f"The answer incorrectly omits the following valid issues: {', '.join(missing_names)}.")
        
        if extra_issues:
            extra_names = [f"'{issues_properties[i]['name']}' (Issue {i})" for i in sorted(list(extra_issues))]
            error_messages.append(f"The answer incorrectly includes the following issues that do not meet the criteria: {', '.join(extra_names)}.")
        
        return f"Incorrect. {' '.join(error_messages)}"

# Execute the check and print the result.
result = check_genomics_question_answer()
print(result)