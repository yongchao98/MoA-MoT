def check_genomics_error_answer():
    """
    Checks the correctness of the answer to a genomics knowledge question.

    The function models the reasoning process by:
    1. Defining each potential issue with properties based on established bioinformatics knowledge.
    2. Identifying which issues meet the question's criteria: "common" AND "difficult-to-spot".
    3. Comparing the derived correct set of issues with the set represented by the given answer.
    """

    # Step 1: Define the knowledge base for each issue.
    # Each issue is evaluated on whether it's a common problem and if the
    # resulting errors are difficult to spot (i.e., they don't cause an obvious crash).
    issues_knowledge_base = {
        1: {
            "name": "Mutually incompatible data formats",
            "is_common": True,
            "is_difficult_to_spot": True,
            "justification": ("While some format errors cause obvious crashes, this category also includes "
                            "classic, subtle errors like 0-based vs. 1-based coordinate confusion. "
                            "This leads to silent off-by-one errors that are very hard to spot.")
        },
        2: {
            "name": "The 'chr' / 'no chr' confusion",
            "is_common": True,
            "is_difficult_to_spot": True,
            "justification": ("A pervasive problem where tools don't crash but silently fail to find "
                            "overlaps, leading to results like '0 variants found' which can be "
                            "misinterpreted as a valid biological finding.")
        },
        3: {
            "name": "Reference assembly mismatch",
            "is_common": True,
            "is_difficult_to_spot": True,
            "justification": ("A fundamental and insidious error. The analysis pipeline runs without "
                            "technical errors, but using data from different genome builds (e.g., hg19 vs hg38) "
                            "makes the results biologically nonsensical. This requires meticulous metadata "
                            "tracking to spot.")
        },
        4: {
            "name": "Incorrect ID conversion",
            "is_common": True,
            "is_difficult_to_spot": True,
            "justification": ("A notorious problem, often caused by spreadsheet software auto-converting "
                            "gene symbols to dates (e.g., SEPT2 -> 'Sep-02'). This silently corrupts data "
                            "before analysis, making it very difficult to trace.")
        }
    }

    # Step 2: Define the options from the question.
    options = {
        "A": {2, 3},
        "B": {1, 2, 3, 4},  # Corresponds to "All of the above"
        "C": {2, 3, 4},
        "D": {3, 4}
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = "B"

    # Step 3: Determine the theoretically correct set of issues based on the knowledge base.
    correct_issues = set()
    for issue_id, properties in issues_knowledge_base.items():
        if properties["is_common"] and properties["is_difficult_to_spot"]:
            correct_issues.add(issue_id)

    # Step 4: Get the set of issues corresponding to the LLM's answer.
    answer_issues = options.get(llm_answer)
    if answer_issues is None:
        return f"Error: The provided answer '{llm_answer}' is not a valid option key."

    # Step 5: Compare the correct set with the answer's set and return the verdict.
    if correct_issues == answer_issues:
        return "Correct"
    else:
        missing_from_answer = correct_issues - answer_issues
        extra_in_answer = answer_issues - correct_issues
        
        error_messages = []
        if missing_from_answer:
            missing_names = [f"'{issues_knowledge_base[i]['name']}'" for i in sorted(list(missing_from_answer))]
            error_messages.append(f"The answer is incorrect because it fails to include the following valid issue(s): {', '.join(missing_names)}.")
        
        if extra_in_answer:
            extra_names = [f"'{issues_knowledge_base[i]['name']}'" for i in sorted(list(extra_in_answer))]
            error_messages.append(f"The answer is incorrect because it includes issue(s) that do not meet the criteria: {', '.join(extra_names)}.")
            
        return "\n".join(error_messages)

# Run the check and print the result.
result = check_genomics_error_answer()
print(result)