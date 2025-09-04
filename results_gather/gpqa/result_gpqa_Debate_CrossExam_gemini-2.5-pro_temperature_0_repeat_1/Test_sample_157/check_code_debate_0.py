def check_biology_question_answer(selected_answer: str):
    """
    Checks the correctness of the answer based on biological principles of dominant-negative mutations.

    Args:
        selected_answer: The letter of the chosen option (e.g., 'A', 'B', 'C', 'D').

    Returns:
        A string indicating "Correct" or the reason for being incorrect.
    """

    # --- Define the core principles from the question ---
    # 1. The protein must form a dimer to be active.
    # 2. The mutation 'Y' is dominant-negative.
    # 3. A dominant-negative mutation in a dimeric protein means the mutant protein
    #    interferes with the wild-type (WT) protein.
    # 4. The most direct way to interfere is for the mutant to still be able to dimerize
    #    with the WT protein, forming a non-functional heterodimer (WT-Y).
    # 5. This interference results in a loss-of-function phenotype, not a gain-of-function.
    # 6. A common cellular fate for misassembled or non-functional protein complexes
    #    is targeted degradation by quality control machinery (e.g., proteasome).
    #    This would destroy both the mutant and the WT protein in the complex.

    # --- Evaluate the provided answer against these principles ---
    if selected_answer == 'A':
        # A) protein aggregation and loss-of-function phenotype
        return "Incorrect: While aggregation can cause a loss-of-function, a more specific and common mechanism for this type of dominant-negative effect is the formation of a non-functional complex that is subsequently targeted for degradation. This actively removes the wild-type protein, which is a more direct explanation for the dominant-negative effect than general aggregation."

    elif selected_answer == 'B':
        # B) change of protein conformation and gain-of-function phenotype
        return "Incorrect: The question explicitly states the mutation is 'dominant-negative', which describes a loss-of-function mechanism, not a gain-of-function."

    elif selected_answer == 'C':
        # C) protein degradation and loss-of-function of the wild-type allele
        # This aligns perfectly with the principles.
        # - The mutant dimerizes with the WT.
        # - The resulting complex is non-functional and recognized as faulty.
        # - The cell's quality control degrades the complex, removing the WT protein.
        # - This causes a loss of function of the WT allele's product.
        return "Correct"

    elif selected_answer == 'D':
        # D) loss of protein dimerization and wild-type phenotype
        return "Incorrect: This option has two flaws. First, if the mutation caused a 'loss of protein dimerization', the mutant protein could not interfere with the wild-type protein, and the effect would not be dominant-negative. Second, the problem states the mutation has a dominant effect, so the phenotype cannot be 'wild-type'."

    else:
        return f"Invalid option '{selected_answer}' provided."

# The answer provided by the LLM is 'C'.
# We run the check on this answer.
result = check_biology_question_answer('C')
print(result)