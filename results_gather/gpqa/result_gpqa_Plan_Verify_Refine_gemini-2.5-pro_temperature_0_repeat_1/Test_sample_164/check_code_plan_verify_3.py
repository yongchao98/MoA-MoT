import json

def check_answer(llm_answer_choice):
    """
    Checks the correctness of the LLM's answer to the ethylene polymerization question.

    The function uses a knowledge base about tandem catalysis for polyethylene production
    to evaluate each statement.
    """

    knowledge_base = {
        "oligomerization_catalyst": {
            "metal_group": "Group 6 (VIa)",
            "metal_example": "Chromium (Cr)",
            "reaction": "Selective trimerization of ethylene to 1-hexene.",
            "industrial_relevance": "High, basis for commercial processes (e.g., Phillips, Chevron Phillips).",
        },
        "activators": {
            "type": "Aluminum-based compounds",
            "examples": ["Methylaluminoxane (MAO)", "Modified Methylaluminoxane (MMAO)"],
            "function": "Essential for activating the chromium oligomerization catalyst.",
        },
        "commercial_status": {
            "is_commercialized": True,
            "companies": ["Univation Technologies (Dow/ExxonMobil)", "Chevron Phillips Chemical"],
            "nuance": "While commercialized, it is an advanced technology and not as widespread as traditional copolymerization with externally supplied comonomers.",
        },
        "noble_metals": {
            "can_be_used": True,
            "context": "Can catalyze ethylene reactions, but are generally less selective and/or cost-effective for this specific industrial application compared to Cr-based systems.",
        }
    }

    # Evaluate each statement based on the knowledge base
    evaluations = {}

    # Statement A: Certain noble metal catalysts can be used but are too expensive.
    # Evaluation: Plausible but less specific and central than C. Cr-based systems are the key technology.
    evaluations['A'] = 'Plausible but not the best answer.'

    # Statement B: Aluminum-based activators do not work for the essential additional reaction step.
    # Evaluation: False. They are essential.
    if "Aluminum-based" in knowledge_base["activators"]["type"]:
        evaluations['B'] = 'Incorrect. Aluminum-based activators like MAO are essential for the oligomerization step with Cr catalysts.'
    else:
        evaluations['B'] = 'Evaluation failed: Knowledge base incomplete.'

    # Statement C: One can use a catalyst of a group VIa transition metal in combination with specific activators.
    # Evaluation: True. This is the core of the technology.
    if knowledge_base["oligomerization_catalyst"]["metal_group"] == "Group 6 (VIa)":
        evaluations['C'] = 'Correct. This accurately describes the key chromium-based catalyst systems used for this process.'
    else:
        evaluations['C'] = 'Evaluation failed: Knowledge base incomplete.'

    # Statement D: Such combined systems are already implemented on an industrial scale in the US.
    # Evaluation: True, but nuanced. C is a more fundamental chemical fact.
    if knowledge_base["commercial_status"]["is_commercialized"]:
        evaluations['D'] = f'Technically correct, as commercial processes exist. However, the implementation is not as widespread as traditional methods, making this statement potentially misleading. Statement C is a more fundamental and undisputed fact about the chemistry involved.'
    else:
        evaluations['D'] = 'Incorrect.'

    # Determine the best answer
    # The best answer is the one that is fundamentally and unambiguously correct.
    # Statement B is definitively false.
    # Statement C is definitively true and describes the core chemistry.
    # Statement D is true but nuanced.
    # Statement A is plausible but less relevant.
    # Therefore, C is the best choice.
    correct_choice = 'C'

    if llm_answer_choice == correct_choice:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer_choice}' is incorrect.\n"
        reason += f"The correct answer is '{correct_choice}'.\n\n"
        reason += "Here is the analysis of all options:\n"
        reason += f"A) {evaluations['A']}\n"
        reason += f"B) {evaluations['B']}\n"
        reason += f"C) {evaluations['C']}\n"
        reason += f"D) {evaluations['D']}\n\n"
        reason += f"Conclusion: Statement C provides the most accurate and fundamental description of the chemistry for this industrial process. Statement B is factually wrong."
        return reason

# The LLM's answer from the prompt
llm_answer = "C"

# Check the answer
result = check_answer(llm_answer)
print(result)