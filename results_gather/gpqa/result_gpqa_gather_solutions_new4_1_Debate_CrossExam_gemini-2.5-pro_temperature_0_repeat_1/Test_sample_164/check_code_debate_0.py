import json

def check_answer_correctness():
    """
    Checks the correctness of the provided LLM answer based on a knowledge base
    of polymer chemistry facts.
    """
    # The final answer from the LLM to be checked.
    llm_answer = "A"

    # A knowledge base representing the ground truth for each statement.
    # This is built from a consensus of chemical literature and industrial practice.
    knowledge_base = {
        "A": {
            "is_correct": True,
            "reason": "This statement is correct and unambiguous. Chromium (Cr), a Group VIa metal, is the key component in state-of-the-art industrial catalysts (e.g., the Chevron Phillips process) for the *selective* trimerization of ethylene to 1-hexene. This specific reaction is essential for creating the 'regular branches' mentioned in the question. These catalysts require specific activators to function."
        },
        "B": {
            "is_correct": False,
            "reason": "This statement is factually incorrect. The 'essential additional reaction step' (selective oligomerization) is catalyzed by systems that almost universally require aluminum-based activators (co-catalysts) like methylaluminoxane (MAO) or trialkylaluminum compounds to function. They are essential, they don't 'not work'."
        },
        "C": {
            "is_correct": False, # Considered incorrect in a single-best-answer context due to ambiguity.
            "reason": "This statement is misleading. While the core technology (selective oligomerization) is industrialized in the US, it's typically in dedicated plants that produce alpha-olefins as a separate product. A true 'one-pot' or single-reactor dual-catalyst system is not the widespread industrial standard. An ambiguous statement is not the best answer in a multiple-choice question."
        },
        "D": {
            "is_correct": False,
            "reason": "This statement is incorrect for the specific goal of the question. Noble metal catalysts (like Palladium) typically produce branched polyethylene via a 'chain-walking' mechanism. This creates a hyperbranched structure with a variety of different branch lengths, not the 'regular branches' of a uniform length that result from incorporating a specific comonomer."
        }
    }

    # Determine the single best answer from the knowledge base.
    # For this problem, there is only one statement that is both correct and unambiguous.
    try:
        correct_answer = [option for option, data in knowledge_base.items() if data["is_correct"]][0]
    except IndexError:
        return "Error: The knowledge base does not contain a correct answer."

    # Compare the LLM's answer to the ground truth.
    if llm_answer == correct_answer:
        return "Correct"
    else:
        llm_reason = knowledge_base.get(llm_answer, {}).get("reason", "The provided option is invalid.")
        correct_reason = knowledge_base.get(correct_answer, {}).get("reason", "")
        
        error_message = (
            f"The provided answer '{llm_answer}' is incorrect.\n"
            f"Reason why '{llm_answer}' is wrong: {llm_reason}\n"
            f"The correct answer is '{correct_answer}'. Reason: {correct_reason}"
        )
        return error_message

# Run the check and print the result.
result = check_answer_correctness()
print(result)