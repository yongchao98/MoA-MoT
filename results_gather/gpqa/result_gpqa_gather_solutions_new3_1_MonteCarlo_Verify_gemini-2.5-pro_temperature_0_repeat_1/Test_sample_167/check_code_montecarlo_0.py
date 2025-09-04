import re

def check_genomics_question_answer(llm_answer_text):
    """
    Checks the correctness of the answer to the genomics data analysis question.

    The function codifies the expert consensus on which issues are both common and
    notoriously difficult to spot (i.e., cause silent failures).
    """

    # --- Knowledge Base ---
    # Define the properties of each issue based on expert consensus in bioinformatics.
    # The critical distinction is "difficult-to-spot," meaning the error doesn't
    # cause an obvious crash but leads to scientifically invalid results.
    issues_properties = {
        1: {
            "name": "Mutually incompatible data formats",
            "is_common": True,
            "is_difficult_to_spot": False,
            "reasoning": "This is a broad category. While some format errors are subtle (e.g., 0-based vs. 1-based coordinates), many common incompatibilities (e.g., wrong file type) cause immediate program crashes, making them easy to spot."
        },
        2: {
            "name": "The 'chr' / 'no chr' confusion",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reasoning": "This is a classic silent failure. The analysis runs but produces zero or incomplete results, which can be misinterpreted as a valid biological finding."
        },
        3: {
            "name": "Reference assembly mismatch",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reasoning": "This is a highly insidious error. The analysis produces plausible-looking but fundamentally incorrect results because coordinates are mapped to the wrong genomic features."
        },
        4: {
            "name": "Incorrect ID conversion",
            "is_common": True,
            "is_difficult_to_spot": True,
            "reasoning": "This error, like Excel auto-converting gene names to dates, silently corrupts the input data. Downstream analysis runs without error but on flawed data."
        }
    }

    # --- Logic ---
    # The question asks for the MOST common sources of DIFFICULT-TO-SPOT errors.
    # We filter our knowledge base for issues that fit this description.
    correct_issue_numbers = {
        num for num, properties in issues_properties.items() if properties["is_difficult_to_spot"]
    }
    # Based on our knowledge base, the correct set is {2, 3, 4}.

    # Map the multiple-choice options to sets of issue numbers.
    option_map = {
        "A": {2, 3},
        "B": {3, 4},
        "C": {1, 2, 3, 4},  # Corresponds to "All of the above"
        "D": {2, 3, 4}
    }

    # --- Verification ---
    # Extract the final answer from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<A>>>, <<<B>>>, etc."
    
    selected_option = match.group(1)
    selected_issues = option_map.get(selected_option)

    if selected_issues == correct_issue_numbers:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrectness.
        missing_issues = correct_issue_numbers - selected_issues
        extra_issues = selected_issues - correct_issue_numbers
        
        reasons = []
        if missing_issues:
            missing_names = [f"'{issues_properties[i]['name']}' (Issue {i})" for i in sorted(list(missing_issues))]
            reasons.append(f"The answer is incorrect because it omits the following common and difficult-to-spot issues: {', '.join(missing_names)}.")
        
        if extra_issues:
            for issue_num in sorted(list(extra_issues)):
                reasoning = issues_properties[issue_num]['reasoning']
                reasons.append(f"The answer is incorrect because it includes '{issues_properties[issue_num]['name']}' (Issue {issue_num}). This issue does not fit the 'difficult-to-spot' criterion as well as others, because: {reasoning}")

        return "\n".join(reasons) if reasons else "Incorrect for an unknown reason."

# The final answer provided by the LLM to be checked.
final_answer_from_llm = """
Here is a step-by-step analysis of the question and the provided candidate answers.

### 1. Deconstruct the Question
The question asks to identify the **most common** sources of **difficult-to-spot erroneous results**. This has two key components:
- **Most Common:** The issue must occur frequently in real-world genomics data analysis.
- **Difficult-to-Spot:** The error must be subtle. It should not cause an obvious program crash or error message. Instead, it should produce results that appear plausible but are scientifically incorrect (a "silent failure").

### 2. Analyze Each Potential Issue
Let's evaluate each of the four listed issues against these two criteria, based on the reasoning provided in the candidate answers.

*   **1. Mutually incompatible data formats:**
    *   **Common?** Yes, all answers agree this is a very common problem.
    *   **Difficult-to-spot?** This is the point of contention.
        *   Some errors in this category are **easy to spot**, such as providing a BAM file to a tool that expects a VCF. The program will almost certainly crash immediately.
        *   However, other errors are **difficult to spot**. A classic example is the 0-based vs. 1-based coordinate system confusion in formats like BED or GFF. A tool might process the file without any error, but all genomic positions will be off by one, leading to incorrect conclusions.
    *   **Conclusion:** This category is a mixed bag. It is a source of both easy-to-spot and difficult-to-spot errors.

*   **2. The "chr" / "no chr" confusion:**
    *   **Common?** Yes, all answers agree this is a pervasive and frequent problem due to lack of standardization between major resources like UCSC and Ensembl.
    *   **Difficult-to-spot?** Yes. This error typically does not cause a crash. Instead, a tool looking for overlaps between two files with different chromosome naming conventions will simply find no matches. It will output "0 results," which a researcher might mistakenly interpret as a true biological negative result rather than a technical failure. This is a canonical example of a silent failure.

*   **3. Reference assembly mismatch:**
    *   **Common?** Yes, this is a very common mistake, especially when combining data from different sources or over time as reference genomes are updated (e.g., hg19 vs. hg38).
    *   **Difficult-to-spot?** Yes, this is one of the most insidious errors. The analysis tools will run to completion, producing results (e.g., variant calls, gene expression levels) that look plausible. However, because the genomic coordinates are being interpreted against the wrong map, the biological conclusions are fundamentally flawed.

*   **4. Incorrect ID conversion:**
    *   **Common?** Yes, this is extremely common and well-documented, with studies showing it affects a large percentage (20-30%) of published papers that use spreadsheet software like Excel.
    *   **Difficult-to-spot?** Yes. The conversion of gene symbols like `SEPT2` to `2-Sep` happens automatically and silently corrupts gene lists. Downstream analyses, like pathway analysis, will run without error but will be based on an incomplete or incorrect set of genes, leading to invalid conclusions.

### 3. Synthesize and Select the Best Answer

We have established that issues 2, 3, and 4 are all common and are almost always difficult to spot, perfectly matching the question's criteria.

The debate centers on issue 1. While it is undoubtedly a source of difficult-to-spot errors (like the coordinate system issue), it is a very broad category that also includes many obvious, easy-to-spot errors (like using the wrong file extension).

The question asks for the "most common sources of *difficult-to-spot* erroneous results." Issues 2, 3, and 4 are more purely representative of this type of silent failure than the mixed-bag category of issue 1. In the context of a multiple-choice question, the option that includes the most precise and characteristic examples is the strongest. Issues 2, 3, and 4 are the "poster children" for insidious, silent errors in bioinformatics.

Therefore, the set {2, 3, 4} is the best and most accurate answer. This corresponds to option D. The majority of the detailed candidate answers (Answers 2, 4, 5, 6, 12) also converge on this conclusion.

<<<D>>>
"""

# Run the check
result = check_genomics_question_answer(final_answer_from_llm)
print(result)