import re

def check_genomics_question_answer(llm_response: str) -> str:
    """
    Checks the correctness of an LLM's answer to a multiple-choice question about genomics data analysis errors.

    The function uses a predefined knowledge base representing expert consensus to determine the correct answer.
    It then checks if the LLM's final choice and its reasoning are consistent with this correct answer.

    Args:
        llm_response: A string containing the LLM's full response, including its reasoning and final answer.

    Returns:
        A string indicating "Correct" or a reason for the answer being incorrect.
    """

    # 1. Define the problem's context and a knowledge base based on expert consensus.
    # The question asks for "common" and "difficult-to-spot" (silent failure) errors.
    # The consensus from detailed analyses is that all four issues can cause such errors.
    knowledge_base = {
        1: {"name": "Mutually incompatible data formats", "is_correct": True, "reason": "Can cause subtle, silent errors like 0-based vs. 1-based coordinate mismatches."},
        2: {"name": "The 'chr' / 'no chr' confusion", "is_correct": True, "reason": "Classic silent failure where tools report 0 overlaps instead of crashing."},
        3: {"name": "Reference assembly mismatch", "is_correct": True, "reason": "Insidious silent failure where analysis runs on mismatched coordinates, producing plausible but wrong results."},
        4: {"name": "Incorrect ID conversion", "is_correct": True, "reason": "Notorious silent failure where gene IDs are corrupted (e.g., by Excel) before analysis begins."}
    }

    # 2. Determine the correct set of issues and the corresponding option letter.
    correct_issue_numbers = {num for num, data in knowledge_base.items() if data["is_correct"]}
    
    options = {
        "A": {3, 4},
        "B": {2, 3},
        "C": {2, 3, 4},
        "D": {1, 2, 3, 4}  # "All of the above"
    }

    correct_option_letter = None
    for letter, issue_set in options.items():
        if issue_set == correct_issue_numbers:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        # This is a sanity check for the checker's own logic.
        return "Checker Error: Could not determine the correct option based on the knowledge base."

    # 3. Extract the LLM's final answer from the text.
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Constraint not satisfied: The final answer was not provided in the required format '<<<...>>>'."
    
    llm_choice = match.group(1)

    # 4. Check if the LLM's choice is correct.
    if llm_choice != correct_option_letter:
        return f"Incorrect. The final answer provided was '{llm_choice}', but the most comprehensive answer is '{correct_option_letter}'. The reasoning should conclude that all four issues are common sources of difficult-to-spot errors."

    # 5. Check if the LLM's reasoning is consistent with its choice.
    # For choice 'D', the reasoning must support that all four issues are valid.
    reasoning_text = llm_response.split("<<<")[0].lower()
    
    # Check if the reasoning for each issue concludes it's a valid source of difficult-to-spot errors.
    if "this category is a valid source" not in reasoning_text:
        return f"Incorrect. The final answer '{llm_choice}' is correct, but the reasoning for issue 1 (data formats) is not conclusive."
    if "quintessential example of a silent failure" not in reasoning_text and "classic and pervasive problem" not in reasoning_text:
        return f"Incorrect. The final answer '{llm_choice}' is correct, but the reasoning for issue 2 ('chr' confusion) is weak."
    if "one of the most insidious errors" not in reasoning_text:
        return f"Incorrect. The final answer '{llm_choice}' is correct, but the reasoning for issue 3 (assembly mismatch) is weak."
    if "notoriously common and widespread problem" not in reasoning_text:
        return f"Incorrect. The final answer '{llm_choice}' is correct, but the reasoning for issue 4 (ID conversion) is weak."

    # Check if the synthesis step correctly concludes that all four are correct.
    if "all four are correct" not in reasoning_text and "all of them are correct" not in reasoning_text:
        return f"Incorrect. The final answer '{llm_choice}' is correct, but the synthesis step does not explicitly state that all four issues are valid, which is the justification for choosing 'All of the above'."

    return "Correct"

# The user-provided response to check
llm_answer_to_check = """
Here is a step-by-step analysis of the question to determine the most accurate answer.

### Step 1: Deconstruct the Question
The question asks to identify the issues that are the **most common sources** of **difficult-to-spot erroneous results** in genomics data analysis. This requires evaluating each option against two key criteria:
1.  **Common**: The issue must be a frequent problem in real-world bioinformatics workflows.
2.  **Difficult-to-spot**: The error should not cause an obvious program crash or a clear error message. Instead, it should lead to a "silent failure," where the analysis completes and produces results that appear plausible but are scientifically incorrect.

### Step 2: Analyze Each Issue Against the Criteria

1.  **Mutually incompatible data formats:**
    *   **Is it common?** Yes, this is an extremely common problem due to the vast number of tools and file formats (BAM, VCF, BED, GTF, etc.) in bioinformatics.
    *   **Is it difficult to spot?** This is the most nuanced point. While a major incompatibility (e.g., providing a text file where a binary file is expected) will cause an immediate, *easy-to-spot* crash, many subtle incompatibilities lead to silent errors. A prime example is the difference between 0-based and 1-based genomic coordinate systems. A tool may process a file with the wrong coordinate system without any error, but all results will be off by one. Similarly, differences in quality score encodings or non-standard fields in a VCF file can be silently misinterpreted. Because these subtle, non-crashing errors are also common, this category is a valid source of difficult-to-spot errors.

2.  **The "chr" / "no chr" confusion:**
    *   **Is it common?** Yes, this is a classic and pervasive problem due to different chromosome naming conventions from major data sources like UCSC ("chr1") and Ensembl ("1").
    *   **Is it difficult to spot?** Yes. This is a quintessential example of a silent failure. A tool looking for overlaps between files with mismatched chromosome names will often run to completion without any error but will simply report "0 overlaps found." A researcher might mistakenly interpret this as a true negative biological result rather than a technical artifact.

3.  **Reference assembly mismatch:**
    *   **Is it common?** Yes, this is a frequent and critical error, especially when combining data generated at different times or from different labs that may have used different reference genome versions (e.g., hg19 vs. hg38).
    *   **Is it difficult to spot?** Yes, this is one of the most insidious errors. The analysis pipeline will run without any warnings, but the genomic coordinates will be mapped to the wrong genes or features, rendering the biological conclusions completely invalid. Detecting this requires meticulous tracking of metadata for all input files.

4.  **Incorrect ID conversion:**
    *   **Is it common?** Yes, this is a notoriously common and widespread problem. The most famous example is spreadsheet software like Microsoft Excel auto-converting gene symbols (e.g., "SEPT2", "MARCH1") into dates ("2-Sep", "1-Mar"). Studies have shown this affects a significant percentage of published papers.
    *   **Is it difficult to spot?** Yes. The error happens during data handling, often before the main analysis. The analysis tools then run on a corrupted gene list, leading to missing data and flawed conclusions (e.g., in a pathway analysis). The error is not in the computation itself, making it very hard to detect.

### Step 3: Synthesize and Conclude
All four issues are undeniably common in genomics data analysis. Crucially, all four can lead to "silent failures" where the analysis completes without crashing but produces scientifically invalid results. Issues 2, 3, and 4 are classic examples of semantic errors where the data's context is wrong. Issue 1 is a broader category that includes both obvious crashes and the same kind of subtle, semantic errors (like coordinate system mismatches). Since all four are well-documented, frequent, and insidious sources of error, the most comprehensive and accurate answer is that all of them are correct.

### Step 4: Select the Final Answer
Based on the analysis that all four listed issues are common sources of difficult-to-spot errors, the correct option is the one that includes all of them.

- A) 3 and 4 - Incomplete.
- B) 2 and 3 - Incomplete.
- C) 2, 3 and 4 - Incomplete, as it omits the valid and common subtle errors from format incompatibilities.
- D) All of the above - This is the most accurate and comprehensive choice.

<<<D>>>
"""

# Run the check
result = check_genomics_question_answer(llm_answer_to_check)
print(result)