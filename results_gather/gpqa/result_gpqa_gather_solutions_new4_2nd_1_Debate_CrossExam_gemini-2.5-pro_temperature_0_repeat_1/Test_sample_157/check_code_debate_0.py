import re

def check_answer_correctness(llm_response):
    """
    Checks the correctness of the LLM's answer to the molecular biology question.

    The function analyzes the core concepts of the question:
    1.  The protein is a transcription factor that must form a dimer to be active.
    2.  The mutation 'Y' is a 'dominant-negative' mutation in the dimerization domain.

    It then evaluates the provided answer based on the implications of these concepts.
    """

    # --- Core Scientific Principles ---

    # 1. Definition of a Dominant-Negative Mutation:
    #   - It's a LOSS-of-function mutation, not gain-of-function.
    #   - It's DOMINANT, meaning it causes a mutant phenotype even with a wild-type (WT) allele present.
    #   - The MECHANISM involves the mutant protein interfering with ("poisoning") the WT protein.

    # 2. Implications for a Dimeric Protein with a mutation in the dimerization domain:
    #   - To interfere, the mutant protein MUST still be able to dimerize with the WT protein.
    #     A complete loss of dimerization would make the mutation recessive, not dominant-negative.
    #   - The resulting mutant-WT heterodimer is non-functional. This sequesters the WT protein.
    #   - This sequestration is a "loss-of-function of the wild-type allele" because the WT protein product is rendered useless.
    #   - A common cellular response to such non-functional or misfolded protein complexes is to target them for degradation via quality control pathways (e.g., the proteasome).

    # --- Evaluation of Options ---

    options_analysis = {
        'A': {
            "is_correct": True,
            "reason": "This is the most accurate and comprehensive description. 'Loss-of-function of the wild-type allele' perfectly describes the genetic effect of a dominant-negative mutation. 'Protein degradation' is a highly plausible and common molecular mechanism for how the cell eliminates the non-functional 'poisoned' complexes."
        },
        'B': {
            "is_correct": False,
            "reason": "This is incorrect because a dominant-negative mutation causes a LOSS-of-function, not a GAIN-of-function."
        },
        'C': {
            "is_correct": False,
            "reason": "This is plausible but less precise than option A. While aggregation can cause a loss of function, it is a specific physical outcome. Option A describes the fundamental genetic effect ('loss-of-function of the wild-type allele') and a more general and direct cellular response (degradation)."
        },
        'D': {
            "is_correct": False,
            "reason": "This is incorrect. For the mutation to be dominant-negative, the mutant protein must be able to dimerize to interfere with the wild-type protein. A complete loss of dimerization would lead to a recessive mutation. Also, a dominant mutation causes a mutant phenotype, not a wild-type one."
        }
    }

    # --- Check the LLM's Answer ---
    
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "Error: Could not find a final answer in the standard format '<<<X>>>' in the provided text."

    llm_choice = match.group(1)

    if options_analysis[llm_choice]["is_correct"]:
        return "Correct"
    else:
        error_message = f"The provided answer <<< {llm_choice} >>> is incorrect.\n"
        error_message += f"Reason: {options_analysis[llm_choice]['reason']}\n\n"
        
        correct_choice = None
        for option, details in options_analysis.items():
            if details["is_correct"]:
                correct_choice = option
                break
        
        error_message += f"The correct answer is <<< {correct_choice} >>>.\n"
        error_message += f"Reason: {options_analysis[correct_choice]['reason']}"
        
        return error_message

# The final answer provided for checking
final_llm_answer = """
Here is a step-by-step analysis of the question to determine the most likely molecular phenotype.

### 1. Deconstruct the Key Concepts
*   **Protein Function:** The protein is a transcription factor that must form a dimer (a complex of two subunits) to be active. This means dimerization is an essential step for its function.
*   **Heterozygous State:** The cell contains one normal (wild-type, WT) allele and one allele with mutation Y. Consequently, the cell produces both normal WT protein subunits and mutant Y protein subunits.
*   **Dominant-Negative Mutation:** This is the most critical concept. A dominant-negative mutation produces a mutant protein that has two key properties:
    *   It is non-functional on its own.
    *   It actively interferes with the function of the normal protein produced from the wild-type allele. This is why the mutation is "dominant" (it has an effect even when a WT allele is present) and "negative" (it causes a loss of function).

### 2. Analyze the Molecular Mechanism of Interference
*   **Location of Mutation:** The mutation (Y) is in the dimerization domain.
*   **Interaction is Required:** For the mutant Y protein to interfere with the WT protein, it must be able to physically interact with it. Since the protein's function and the mutation's location are related to dimerization, the interference must occur at this step.
*   **Formation of "Poisoned" Dimers:** The mutant Y subunit must still be able to bind to a WT subunit. This forms a mixed dimer, or heterodimer (WT-Y). If the mutation completely prevented dimerization, the Y protein would be inert and could not interfere with the WT proteins, which would be free to form functional WT-WT dimers. This would result in a recessive mutation, not a dominant-negative one.
*   **Consequence:** The resulting WT-Y heterodimer is non-functional. The mutant Y subunit effectively "poisons" the complex, trapping the functional WT subunit in an inactive state. This sequestration of functional protein is the core of the dominant-negative effect.

### 3. Evaluate the Options Based on the Mechanism
*   **B) change of protein conformation and gain-of-function phenotype:** This is incorrect. A dominant-*negative* mutation is a type of loss-of-function mutation, not a gain-of-function.
*   **D) loss of protein dimerization and wild-type phenotype:** This is incorrect. As explained above, the mutant protein *must* be able to dimerize to interfere with the wild-type protein. A complete loss of dimerization would lead to a recessive mutation. Furthermore, a dominant mutation causes a mutant phenotype, not a wild-type one.
*   **C) protein aggregation and loss-of-function phenotype:** This is a plausible but less precise answer. The outcome is indeed a loss-of-function. A mutation in an interaction domain could lead to misfolding and cause proteins to form large, non-specific clumps (aggregates), which would trap the WT protein. However, the primary mechanism is the formation of a specific, non-functional *dimer*. Aggregation is a possible secondary consequence, but not the most fundamental description.
*   **A) protein degradation and loss-of-function of the wild-type allele:** This is the most accurate and comprehensive description.
    *   **"loss-of-function of the wild-type allele"**: This phrase perfectly describes the genetic outcome of a dominant-negative mutation. The protein produced by the normal allele is rendered non-functional by the mutant protein's interference.
    *   **"protein degradation"**: This is a very common and well-documented molecular fate for aberrant or non-functional protein complexes. The cell's quality control machinery (like the ubiquitin-proteasome system) recognizes the faulty WT-Y heterodimer and targets the entire complex for destruction. This degradation is an active mechanism that removes the "poisoned" wild-type protein from the cell, cementing the loss of function.

### Conclusion
The defining characteristic of the dominant-negative mutation Y is its ability to form a non-functional complex with the wild-type protein, thereby inhibiting its function. This is precisely described as a "loss-of-function of the wild-type allele." A highly probable cellular response to such faulty complexes is to target them for "protein degradation." Therefore, this option provides the most complete and accurate description of the molecular phenotype.

<<<A>>>
"""

# Run the check
result = check_answer_correctness(final_llm_answer)
print(result)