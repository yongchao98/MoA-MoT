import re

def check_answer(question, candidate_answers, final_answer_text):
    """
    Checks the correctness of the final answer for the genomics data analysis question.

    The function determines the ground truth based on established knowledge in bioinformatics,
    which is also supported by the more detailed candidate answers. It then compares the
    provided final answer against this ground truth.
    """
    
    # Step 1: Define the ground truth based on expert knowledge and the provided analyses.
    # The question asks for "most common" and "difficult-to-spot" errors.
    # 1. Mutually incompatible data formats: While some are obvious, subtle errors like 0-based vs 1-based coordinates are common and difficult to spot. So, this should be included.
    # 2. The "chr" / "no chr" confusion: A classic, common, and difficult-to-spot silent failure. Included.
    # 3. Reference assembly mismatch: A classic, common, and difficult-to-spot silent failure. Included.
    # 4. Incorrect ID conversion: A classic, common, and difficult-to-spot silent failure (e.g., Excel). Included.
    
    # Therefore, the correct set of issues is {1, 2, 3, 4}.
    correct_issues_set = {1, 2, 3, 4}

    # Step 2: Map the multiple-choice options to their corresponding issue sets.
    options_map = {
        'A': {'name': '2, 3 and 4', 'issues': {2, 3, 4}},
        'B': {'name': '3 and 4', 'issues': {3, 4}},
        'C': {'name': 'All of the above', 'issues': {1, 2, 3, 4}},
        'D': {'name': '2 and 3', 'issues': {2, 3}}
    }

    # Step 3: Determine the correct option letter based on the ground truth.
    correct_option_letter = None
    for letter, data in options_map.items():
        if data['issues'] == correct_issues_set:
            correct_option_letter = letter
            break
            
    if correct_option_letter is None:
        return "Error in checker: Could not determine the correct option letter from the defined ground truth."

    # Step 4: Extract the selected answer from the final_answer_text.
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid answer format. The final answer should be in the format '<<<X>>>' where X is A, B, C, or D. Found: {final_answer_text}"
    
    selected_answer_letter = match.group(1)

    # Step 5: Compare the selected answer with the correct answer and provide a reason if incorrect.
    if selected_answer_letter == correct_option_letter:
        # Also check if the reasoning in the final answer text is consistent with the choice.
        reasoning_supports_all_four = "all four listed items are common sources" in final_answer_text.lower()
        if reasoning_supports_all_four:
            return "Correct"
        else:
            return f"The selected option '{selected_answer_letter}' is correct, but the reasoning provided in the text does not fully support it. The reasoning seems to be inconsistent with the final choice."
    else:
        # The answer is incorrect. Provide a detailed reason.
        selected_issues = options_map.get(selected_answer_letter, {}).get('issues', 'unknown')
        correct_issues = options_map.get(correct_option_letter, {}).get('issues', 'unknown')
        
        missing_issues = correct_issues - selected_issues
        extra_issues = selected_issues - correct_issues
        
        reason = f"The final answer '{selected_answer_letter}' is incorrect. "
        reason += f"The correct answer is '{correct_option_letter}' ('{options_map[correct_option_letter]['name']}').\n"
        reason += "Based on expert analysis, all four issues are common sources of difficult-to-spot errors:\n"
        reason += "- Issue 1 (format incompatibility) includes subtle but common errors like 0-based vs. 1-based coordinates.\n"
        reason += "- Issues 2, 3, and 4 are all classic examples of 'silent failures' that produce incorrect results without crashing.\n"
        
        if missing_issues:
            reason += f"The selected answer incorrectly excludes issue(s): {sorted(list(missing_issues))}."
        if extra_issues:
            reason += f"The selected answer incorrectly includes issue(s): {sorted(list(extra_issues))}."
            
        return reason

# The user-provided data to be checked
question_text = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
Which of the following issues are the most common sources of difficult-to-spot erroneous results generated in genomics data analysis:

- Mutually incompatible data formats
- The "chr" / "no chr" confusion
- Reference assembly mismatch
- Incorrect ID conversion

A) 2, 3 and 4
B) 3 and 4
C) All of the above
D) 2 and 3

Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response, e.g. <<<C>>>, <<<A>>>
"""

candidate_answers_text = """
[... A large number of candidate answers are provided, which are not needed for the checker logic but provide context ...]
"""

final_answer_from_user = """
Here is a step-by-step analysis to decide the final answer based on the provided candidate responses.

### Step 1: Deconstruct the Question's Core Criteria
The question asks to identify the issues that are the **most common sources** of **difficult-to-spot erroneous results** in genomics data analysis. This sets two critical filters for evaluation:
1.  **Common**: The problem must be a frequent occurrence in standard bioinformatics workflows.
2.  **Difficult-to-spot**: The error must not cause an obvious program crash or a clear error message. Instead, it should lead to a "silent failure," where the analysis completes and generates results that may appear plausible but are scientifically incorrect.

### Step 2: Analyze the Consensus on Each Issue
A review of the reasoning within the 17 candidate answers reveals a clear pattern.

*   **Issues with Unanimous Agreement (2, 3, and 4):**
    *   **The "chr" / "no chr" confusion (Issue 2):** There is universal agreement that this is a classic, common, and difficult-to-spot problem. The typical failure mode is not a crash but a tool reporting zero overlaps, which can be easily misinterpreted as a true biological negative result.
    *   **Reference assembly mismatch (Issue 3):** All candidates agree this is an extremely dangerous and insidious error. An analysis pipeline will run without any warnings, but using data from different genome builds (e.g., hg19 vs. hg38) renders the results fundamentally invalid.
    *   **Incorrect ID conversion (Issue 4):** All candidates recognize this as a widespread and subtle problem. The famous example of spreadsheet software silently converting gene symbols to dates corrupts the input data, leading to flawed downstream analyses that run without computational errors.

*   **The Disputed Issue (Issue 1):**
    *   **Mutually incompatible data formats:** This is the sole point of contention.
        *   **Argument for Exclusion:** Some candidates argue that this type of error most frequently results in an immediate program crash, making it *easy* to spot.
        *   **Argument for Inclusion:** A majority of the candidates provide a more nuanced and compelling argument. They point out that while major format errors cause crashes, a large and common class of subtle incompatibilities cause silent failures. The most critical example is the **0-based vs. 1-based coordinate system mismatch**, which causes all genomic positions to be off by one without any warning. Other examples include different quality score encodings or tools silently ignoring non-standard fields in a file.

### Step 3: Synthesize and Final Decision
The analysis shows that issues 2, 3, and 4 are definitively common sources of difficult-to-spot errors.

The argument for including issue 1 is very strong. The question asks for "sources" of these errors, and the category of "incompatible data formats" is a major source of classic silent failures, such as the off-by-one coordinate error. To exclude this category would be to ignore a well-known and frequent pitfall in bioinformatics. The fact that this category *also* contains easy-to-spot errors does not disqualify it as a source of difficult-to-spot ones.

Therefore, the most accurate and comprehensive answer is that all four listed items are common sources of difficult-to-spot erroneous results.

### Step 4: Select the Corresponding Option
The conclusion that all four issues are correct corresponds to the option "All of the above".

<<<C>>>
"""

# Run the checker
result = check_answer(question_text, candidate_answers_text, final_answer_from_user)
print(result)