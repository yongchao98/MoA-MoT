import re

def check_answer_correctness(question_details, final_answer):
    """
    Checks the correctness of the answer to a conceptual biology question.

    Args:
        question_details (dict): A dictionary containing the premises of the question.
        final_answer (str): The letter of the chosen answer (e.g., "B").

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the error.
    """

    # --- Step 1: Define the logical constraints from the question ---

    # Constraint 1: From "recessive loss-of-function mutation X", we deduce the protein is haplosufficient.
    # This means 50% of the protein function (from one good allele) is enough for a normal phenotype.
    is_haplosufficient = True

    # Constraint 2: From "dominant-negative mutation Y", we know two things:
    # 2a. The phenotype must be a "loss-of-function".
    # 2b. The phenotype must be "dominant" (i.e., not wild-type in a heterozygote).
    # 2c. The mechanism must involve the mutant protein interfering with the wild-type protein.
    #    A simple loss of the mutant protein's function is not enough to cause a dominant phenotype
    #    due to haplosufficiency.

    # Constraint 3: From "mutation Y in the dimerization domain", the interference mechanism
    # must involve protein-protein interaction (dimerization). The mutant protein must still be
    # able to interact with the wild-type protein to interfere with it.

    # --- Step 2: Define the options provided in the context ---
    # Based on the provided text (e.g., Answer 3), the options are:
    options = {
        "A": "protein aggregation and loss-of-function phenotype",
        "B": "protein degradation and loss-of-function of the wild-type allele",
        "C": "loss of protein dimerization and wild-type phenotype",
        "D": "change of protein conformation and gain-of-function phenotype"
    }

    if final_answer not in options:
        return f"Invalid answer format. The answer should be one of {list(options.keys())}."

    selected_option_text = options[final_answer]

    # --- Step 3: Evaluate the selected answer against the constraints ---

    # Check against Constraint 2a (must be loss-of-function)
    if "gain-of-function" in selected_option_text:
        return "Incorrect. The mutation is defined as dominant-negative, which is a type of loss-of-function, not a gain-of-function."

    # Check against Constraint 2b (must be dominant, not wild-type)
    if "wild-type phenotype" in selected_option_text:
        return "Incorrect. The mutation is defined as dominant, meaning it produces a mutant phenotype in a heterozygote, not a wild-type phenotype."

    # Check against Constraint 3 (must allow for interference)
    if "loss of protein dimerization" in selected_option_text:
        return "Incorrect. If the mutation caused a complete loss of dimerization, the mutant protein could not interact with and interfere with the wild-type protein. Due to haplosufficiency, this would result in a wild-type phenotype, making the mutation recessive, not dominant."

    # At this point, only options A and B remain. Let's analyze them.
    # Both describe a loss-of-function phenotype. We need to find the *most likely* molecular phenotype.

    # Option A: "protein aggregation and loss-of-function phenotype"
    # This is a plausible mechanism for loss-of-function.

    # Option B: "protein degradation and loss-of-function of the wild-type allele"
    # This option is more precise and superior for two reasons:
    # 1. "loss-of-function of the wild-type allele" is the textbook definition of a dominant-negative effect.
    # 2. "protein degradation" is a very common and specific cellular quality control mechanism for non-functional protein complexes (like the poisoned wild-type/mutant dimer).

    if final_answer == "B":
        return "Correct"
    else:
        return f"Incorrect. The chosen answer '{final_answer}' is not the best description. While '{selected_option_text}' might be a possible outcome, option B ('protein degradation and loss-of-function of the wild-type allele') is more precise. It correctly identifies the genetic effect (functional loss of the wild-type allele's product) and proposes a highly plausible molecular mechanism (degradation of the non-functional complex by cellular quality control)."


# --- Execution ---
# The final answer from the LLM analysis is "B".
llm_final_answer = "B"

# The question details can be summarized for the function.
question_info = {
    "protein_type": "dimeric_transcription_factor",
    "mutation_X": {"type": "recessive_loss_of_function"},
    "mutation_Y": {"type": "dominant_negative", "location": "dimerization_domain"}
}

# Run the check.
result = check_answer_correctness(question_info, llm_final_answer)
print(result)