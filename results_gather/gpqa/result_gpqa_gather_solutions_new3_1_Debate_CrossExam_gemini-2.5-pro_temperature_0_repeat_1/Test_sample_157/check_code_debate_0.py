def check_answer_correctness(llm_answer: str):
    """
    Checks the correctness of the LLM's answer based on the principles of molecular genetics.

    Args:
        llm_answer: The answer provided by the LLM (e.g., "A", "B", "C", "D").

    Returns:
        A string indicating "Correct" or the reason for being incorrect.
    """
    # Define the options as described in the problem.
    options = {
        "A": "protein aggregation and loss-of-function phenotype",
        "B": "protein degradation and loss-of-function of the wild-type allele",
        "C": "loss of protein dimerization and wild-type phenotype",
        "D": "change of protein conformation and gain-of-function phenotype"
    }

    if llm_answer not in options:
        return f"Invalid answer format. The answer should be one of {list(options.keys())}."

    selected_text = options[llm_answer]

    # --- Rule 1: Check for contradiction with the "dominant" aspect ---
    # A dominant mutation must produce a mutant phenotype, not a wild-type one.
    if "wild-type phenotype" in selected_text:
        return (f"Incorrect. The answer '{llm_answer}' suggests a 'wild-type phenotype'. "
                f"However, the question states the mutation is dominant, which by definition "
                f"causes a mutant phenotype in a heterozygote.")

    # --- Rule 2: Check for contradiction with the "negative" (loss-of-function) aspect ---
    # A dominant-negative mutation is a type of loss-of-function mutation.
    if "gain-of-function" in selected_text:
        return (f"Incorrect. The answer '{llm_answer}' suggests a 'gain-of-function phenotype'. "
                f"A dominant-negative mutation is a specific class of loss-of-function mutation.")

    # --- Rule 3: Check for contradiction with the dominant-negative mechanism in a dimeric protein ---
    # The mutant protein must dimerize with the wild-type to interfere.
    if "loss of protein dimerization" in selected_text:
        return (f"Incorrect. The answer '{llm_answer}' suggests a 'loss of protein dimerization'. "
                f"For a dominant-negative effect to occur in a dimeric protein, the mutant subunit "
                f"must be able to dimerize with the wild-type subunit to interfere with it. "
                f"A complete loss of dimerization would lead to a recessive mutation.")

    # --- Rule 4: Evaluate the most accurate description ---
    # At this point, only A and B are left. We need to determine which is better.
    # Option B contains two highly accurate concepts:
    # 1. "loss-of-function of the wild-type allele": This is the precise definition of a dominant-negative effect.
    # 2. "protein degradation": This is a common and specific molecular mechanism for how cells handle
    #    the non-functional heterodimers, thus causing the loss of the wild-type protein.
    
    if llm_answer == "B":
        return "Correct"
    
    if llm_answer == "A":
        return ("Incorrect. While option A correctly identifies a 'loss-of-function phenotype', "
                "option B provides a more precise and complete description. 'Loss-of-function of the wild-type allele' "
                "is the definition of the dominant-negative effect, and 'protein degradation' is a common molecular "
                "mechanism for this process. 'Protein aggregation' is a possible but less specific mechanism.")

    # This part should not be reached if the input is A, B, C, or D.
    return "Error in checking logic."

# The final answer from the LLM analysis is 'B'.
# Let's run the check on this answer.
final_answer_to_check = "B"
result = check_answer_correctness(final_answer_to_check)
print(result)