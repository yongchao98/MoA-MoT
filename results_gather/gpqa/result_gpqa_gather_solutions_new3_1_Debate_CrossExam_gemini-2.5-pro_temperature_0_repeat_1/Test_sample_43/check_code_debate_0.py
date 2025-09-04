def check_answer_correctness():
    """
    Checks the correctness of the answer to the mitochondrial experiment question.

    The function encodes the scientific rationale for each option and compares it
    against the provided answer.
    """

    # The final answer provided by the LLM analysis.
    llm_answer = "A"

    # Scientific ground truth for each option.
    # The question asks which experiment will NOT help.
    experimental_validity = {
        "A": {
            "is_helpful": False,
            "reason": "This experiment is NOT helpful. Mitochondria do not directly import glucose; they import pyruvate, the end-product of glycolysis which occurs in the cytoplasm. Therefore, a glucose uptake assay on isolated mitochondria is biochemically incorrect and will not yield useful data."
        },
        "B": {
            "is_helpful": True,
            "reason": "This experiment IS helpful. The chemical name is for JC-1, a standard dye used to measure mitochondrial membrane potential, which is a key indicator of mitochondrial health and function."
        },
        "C": {
            "is_helpful": True,
            "reason": "This experiment IS helpful. It measures extracellular ATP. While an indirect measure, a drug causing severe mitochondrial damage could lead to cell death and ATP release, making this a valid, albeit indirect, method to screen for mitochondrial toxicity."
        },
        "D": {
            "is_helpful": True,
            "reason": "This experiment IS helpful. Mito-RTP is a fluorescent probe that measures mitochondrial temperature. Since mitochondrial respiration generates heat, this provides a direct assessment of metabolic activity."
        }
    }

    # The correct answer is the option that is NOT helpful.
    correct_key = None
    for key, value in experimental_validity.items():
        if not value["is_helpful"]:
            correct_key = key
            break

    if llm_answer == correct_key:
        return "Correct"
    else:
        # The LLM chose an answer that is actually a helpful experiment.
        reason = f"Incorrect. The provided answer '{llm_answer}' is wrong because the experiment described in option {llm_answer} would actually be helpful for investigating the drug's effect on mitochondria. "
        reason += experimental_validity[llm_answer]["reason"]
        return reason

# Execute the check and print the result.
result = check_answer_correctness()
print(result)