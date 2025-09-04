import re

def check_biology_question_answer(question: str, llm_answer: str) -> str:
    """
    Checks the correctness of an answer to a molecular biology question about mutations.

    The function codifies the key principles of dominant-negative mutations in dimeric proteins
    and evaluates the provided answer against these principles.

    Args:
        question: The question text.
        llm_answer: The full response from the LLM, including reasoning and the final answer in <<<X>>> format.

    Returns:
        "Correct" if the answer is correct and well-reasoned.
        A string explaining the error if the answer is incorrect.
    """

    # --- Define the core biological principles for this specific problem ---

    # Principle 1: A dominant-negative mutation is a type of loss-of-function.
    def check_loss_of_function(option_text):
        return "gain-of-function" not in option_text

    # Principle 2: For a dominant-negative effect in a dimeric protein, the mutant must interact with the wild-type.
    # A complete loss of dimerization would lead to a recessive phenotype.
    def check_interaction_is_required(option_text):
        return "loss of protein dimerization" not in option_text

    # Principle 3: The most accurate description includes both the genetic effect and a plausible molecular mechanism.
    # Genetic effect: "loss-of-function of the wild-type allele"
    # Plausible mechanism: "protein degradation" of the faulty complex.
    def check_best_description(option_text):
        return "protein degradation" in option_text and "loss-of-function of the wild-type allele" in option_text

    # --- Extract the chosen answer from the LLM's response ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Error: The answer format is invalid. No answer found in <<<X>>> format."
    
    chosen_option_letter = match.group(1)

    # --- Map option letters to their text content from the question prompt ---
    # This is a simplified representation for the check.
    options = {
        'A': "protein aggregation and loss-of-function phenotype",
        'B': "change of protein conformation and gain-of-function phenotype",
        'C': "loss of protein dimerization and wild-type phenotype",
        'D': "protein degradation and loss-of-function of the wild-type allele"
    }
    
    chosen_option_text = options.get(chosen_option_letter)
    if not chosen_option_text:
        return f"Error: Invalid option '{chosen_option_letter}' selected."

    # --- Evaluate the chosen option against the principles ---

    # Check against Principle 1 (Must be Loss-of-Function)
    if not check_loss_of_function(chosen_option_text):
        return (f"Incorrect. The chosen answer '{chosen_option_letter}' suggests a 'gain-of-function' phenotype. "
                "A dominant-negative mutation is a type of loss-of-function mutation, as it interferes with and reduces the overall protein activity.")

    # Check against Principle 2 (Interaction is Required)
    if not check_interaction_is_required(chosen_option_text):
        return (f"Incorrect. The chosen answer '{chosen_option_letter}' suggests a 'loss of protein dimerization'. "
                "If the mutant protein could not dimerize, it could not interfere with the wild-type protein. "
                "This would result in a recessive mutation, not a dominant-negative one.")

    # Check against Principle 3 (Is it the best, most complete description?)
    if not check_best_description(chosen_option_text):
        # This handles option A, which is plausible but less precise than D.
        if chosen_option_letter == 'A':
            return (f"Incorrect. The chosen answer is '{chosen_option_letter}'. While 'loss-of-function' is correct, 'protein aggregation' is a less specific mechanism. "
                    "Option D is superior because it describes both the precise genetic effect ('loss-of-function of the wild-type allele') "
                    "and a very common molecular mechanism ('protein degradation') for how cells handle such faulty protein complexes.")
        else:
            # This case should not be reached if logic is sound, but as a fallback:
            return f"Incorrect. The chosen answer '{chosen_option_letter}' is not the most accurate description. The correct answer is D."

    # If all checks pass, the answer must be D.
    if chosen_option_letter == 'D':
        return "Correct"
    else:
        # Fallback for any unhandled logic paths.
        return f"An unexpected error occurred in the evaluation logic for option {chosen_option_letter}."


# --- Provided data from the user ---
question_text = """
In an inactive state, a transcription factor subunit awaits an activation signal from the membrane phosphorylation cascade. Upon phosphorylation on the Ser residues in the transactivation domain, it dimerizes, shuttles into the nucleus, and serves in gene transcription. A missense mutation X in the transactivation domain is a recessive loss-of-function mutation. You observe that a heterozygous mutation Y in the dimerization domain acts as a dominant-negative mutation. Which of the following molecular phenotypes is most likely observed in the presence of mutation Y?

A) protein aggregation and loss-of-function phenotype
B) change of protein conformation and gain-of-function phenotype
C) loss of protein dimerization and wild-type phenotype
D) protein degradation and loss-of-function of the wild-type allele
"""

llm_answer_text = """
Here is a step-by-step analysis of the problem to determine the best answer.

**1. Deconstruct the Key Terms**

*   **Transcription Factor:** A protein that controls the rate of transcription of genetic information from DNA to messenger RNA.
*   **Dimerizes:** The protein must pair up with another subunit to become functional. This is a critical step.
*   **Heterozygous Mutation Y:** The cell has one normal (wild-type) allele and one mutated (Y) allele. Therefore, it produces both normal protein subunits and mutant Y protein subunits.
*   **Dimerization Domain:** The mutation is in the part of the protein responsible for pairing up.
*   **Dominant-Negative Mutation:** This is the most important concept. It means the mutant protein is not just non-functional itself, but it also actively interferes with the function of the normal protein produced from the wild-type allele. This results in a loss-of-function phenotype even when a good copy of the gene is present.

**2. Synthesize the Molecular Mechanism**

1.  In a heterozygous cell, both wild-type (WT) and mutant (Y) protein subunits are present.
2.  For the mutation to be "dominant-negative" in a dimeric protein, the mutant Y subunit must still be able to interact with the WT subunit. If it couldn't interact at all, it wouldn't be able to interfere, and the WT subunits could form functional WT-WT dimers, likely resulting in a normal or near-normal phenotype (a recessive mutation).
3.  Therefore, the Y subunit must be able to form a dimer with a WT subunit, creating a WT-Y "heterodimer".
4.  This WT-Y heterodimer is non-functional. The mutant Y subunit acts as a "poison pill," inactivating the entire complex.
5.  By trapping functional WT subunits in these useless complexes, the Y protein drastically reduces the overall activity of the transcription factor. This interference is the essence of the dominant-negative effect.

**3. Evaluate the Answer Choices**

*   **A) protein aggregation and loss-of-function phenotype:** The "loss-of-function phenotype" is correct. However, "protein aggregation" usually implies the formation of large, non-specific, insoluble clumps. While this can happen with misfolded proteins, the more precise mechanism here is the formation of specific, non-functional *dimers*. This option is plausible but less specific than option D.

*   **B) change of protein conformation and gain-of-function phenotype:** This is incorrect. A dominant-negative mutation is a type of *loss-of-function* mutation, not a gain-of-function.

*   **C) loss of protein dimerization and wild-type phenotype:** This is incorrect. If the mutation caused a complete loss of dimerization ability, the mutant protein could not interfere with the wild-type protein. This would lead to a *recessive* mutation, not a dominant-negative one.

*   **D) protein degradation and loss-of-function of the wild-type allele:** This is the most accurate and complete answer.
    *   **"loss-of-function of the wild-type allele"**: This is the very definition of a dominant-negative effect. The mutant protein's presence causes the protein product from the normal allele to become non-functional.
    *   **"protein degradation"**: This is a highly plausible and common molecular mechanism for this effect. Cellular quality control systems (like the ubiquitin-proteasome system) often recognize aberrant or non-functional protein complexes (like the WT-Y heterodimer) and target them for destruction. By causing the degradation of the WT protein it is bound to, the mutant protein effectively eliminates the functional protein from the cell.

**Conclusion**

The presence of a dominant-negative mutation in the dimerization domain leads to the formation of non-functional heterodimers. A common cellular response is to degrade these faulty complexes, which removes the functional wild-type protein along with the mutant one. This mechanism perfectly explains the "loss-of-function of the wild-type allele."

<<<D>>>
"""

# Run the check
result = check_biology_question_answer(question_text, llm_answer_text)
print(result)