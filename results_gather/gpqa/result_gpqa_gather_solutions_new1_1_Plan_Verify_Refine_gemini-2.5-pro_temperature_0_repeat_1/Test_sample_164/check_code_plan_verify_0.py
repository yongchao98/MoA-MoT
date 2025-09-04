import re

def check_correctness_of_answer(question_text, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the ethylene polymerization question.

    The function encodes established chemical facts as a knowledge base and evaluates
    each statement from the original question to determine the single best answer.
    It then compares this to the final answer provided by the LLM.

    Args:
        question_text (str): The original question text.
        llm_answer_text (str): The full text of the LLM's response, including its reasoning and final answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Knowledge Base ---
    # These facts are based on established industrial chemistry for polyolefins.
    knowledge_base = {
        "oligomerization_catalyst_metal": "Chromium",  # Key metal for selective oligomerization
        "metal_group": {
            "Chromium": "VIa"  # or 6
        },
        "oligomerization_activators": ["MAO", "MMAO", "Trialkylaluminum"],  # Common activators
        "activator_base_element": "Aluminum",
        "industrial_implementation": {
            "oligomerization_process": True,  # e.g., Chevron Phillips 1-hexene process
            "location": "US",
            "one_pot_tandem_system": "Rare/Debatable"  # Industrial standard is often separate processes
        },
        "noble_metal_catalysis": {
            "can_be_used": True,
            "mechanism": "chain-walking",  # Different from comonomer incorporation
            "branch_type": "hyperbranched",  # Not "regular"
            "cost": "high",
            "industrial_choice": False
        }
    }

    # --- Evaluation Logic ---
    # We evaluate each statement from the original question to find the most accurate one.
    evaluation = {}
    reasons = {}

    # Evaluate Statement A: "Such combined systems are already implemented on an industrial scale in the US."
    # The term "combined systems" is ambiguous. While the core technology (selective oligomerization)
    # is industrialized in the US, true one-pot tandem systems are not the widespread industrial standard.
    # This makes the statement's correctness debatable and not the best answer.
    evaluation["A"] = "Ambiguous"
    reasons["A"] = "The term 'combined systems' is ambiguous. A strict interpretation (one-pot reactor) is not widely implemented industrially, making the statement debatable."

    # Evaluate Statement B: "Certain noble metal catalysts can be used but are too expensive."
    # While technically true, this is a distracting statement. Noble metals are not the primary choice
    # for creating 'regular branches' via comonomer incorporation.
    evaluation["B"] = "Weak/Distracting"
    reasons["B"] = "This is not the primary technology for this specific application. Noble metals typically produce different polymer architectures and are not the industrial choice for creating regular branches."

    # Evaluate Statement C: "One can use a catalyst of a group VIa transition metal in combination with specific activators."
    # This is a fundamentally correct and central statement. Chromium (Cr), a Group VIa metal, is the key
    # component in state-of-the-art industrial catalysts for the selective oligomerization of ethylene.
    evaluation["C"] = "Correct"
    reasons["C"] = "This is a direct and accurate chemical fact. Chromium (Group VIa) catalysts are the cornerstone of the technology for the 'essential additional reaction step' (selective oligomerization)."

    # Evaluate Statement D: "Aluminum-based activators do not work for the essential additional reaction step."
    # This statement is factually incorrect. The oligomerization catalysts almost universally require
    # an aluminum-based activator (co-catalyst) such as MAO or other organoaluminum compounds.
    evaluation["D"] = "Incorrect"
    reasons["D"] = "This is factually incorrect. Aluminum-based activators (like MAO) are essential for the function of the oligomerization catalysts (e.g., Chromium-based systems)."

    # --- Determine the Best Answer based on our evaluation ---
    correct_answers = [k for k, v in evaluation.items() if v == "Correct"]
    
    if len(correct_answers) == 1:
        best_answer = correct_answers[0]
    else:
        # This case would indicate a flawed question or evaluation logic.
        return f"Error in evaluation logic: Found {len(correct_answers)} unambiguously correct answers. Evaluation: {evaluation}"

    # --- Check the LLM's Provided Answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<X>>> in the provided text."
    
    llm_final_answer = match.group(1)

    if llm_final_answer == best_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is <<<{llm_final_answer}>>>, but the correct answer based on chemical facts is <<<{best_answer}>>>.\n"
                f"Reasoning:\n"
                f"- The LLM chose statement {llm_final_answer}, which is evaluated as '{evaluation[llm_final_answer]}'. Reason: {reasons[llm_final_answer]}\n"
                f"- The correct statement is {best_answer}, which is evaluated as '{evaluation[best_answer]}'. Reason: {reasons[best_answer]}")

# --- Input Data ---
question = """
In an industrial research lab, a scientist performs ethylene polymerization with a homogeneous organometallic catalyst system, generating a polymer of high density. He intends to add a second catalyst system to introduce regular branches in the polymer backbone, also only using ethylene as the reactant. He consults a senior scientist, who gives the following statements. “Such combined systems are already implemented on an industrial scale in the US. One can use a catalyst of a group VIa transition metal in combination with specific activators. Aluminum-based activators do not work for the essential additional reaction step. Certain noble metal catalysts can be used but are too expensive.”
Which of these four statements is correct regarding the formation of a polymer with regular branches using only ethylene as the monomer and a dual catalyst system?

A) Such combined systems are already implemented on an industrial scale in the US.
B) Certain noble metal catalysts can be used but are too expensive.
C) One can use a catalyst of a group VIa transition metal in combination with specific activators.
D) Aluminum-based activators do not work for the essential additional reaction step.
"""

llm_answer = """
Here is a step-by-step analysis to determine the correct answer based on the provided candidate responses.

**Step 1: Deconstruct the Core Chemical Problem**
The task is to create a polymer with "regular branches" using only ethylene as the reactant in a dual-catalyst system. This process is known as tandem catalysis. It requires two distinct reactions to occur:
1.  **Oligomerization:** One catalyst must selectively convert some of the ethylene into a specific short-chain alpha-olefin (e.g., trimerizing three ethylene molecules into one 1-hexene molecule). This is the "essential additional reaction step."
2.  **Copolymerization:** A second catalyst must then polymerize the remaining ethylene with the newly formed alpha-olefin to create a branched polymer chain (Linear Low-Density Polyethylene, or LLDPE). The "regular branches" result from the high selectivity of the oligomerization catalyst.

**Step 2: Analyze the Consensus Reasoning for Each Statement**
Despite significant confusion in the lettering of the final answers across the candidates, the underlying reasoning shows a strong consensus on the chemical facts.

*   **Statement A: "Such combined systems are already implemented on an industrial scale in the US."**
    *   **Analysis:** The candidate answers are divided on this point. While the core technologies (like producing 1-hexene from ethylene) are indeed industrialized in the US (e.g., by Chevron Phillips, Sasol), the term "combined system" is ambiguous. If it implies a single reactor where both oligomerization and polymerization occur simultaneously, its large-scale industrial implementation is not as widespread as the separate production of the components. Therefore, the accuracy of this statement is debatable and depends on interpretation.

*   **Statement B: "Certain noble metal catalysts can be used but are too expensive."**
    *   **Analysis:** The consensus is that while this statement is generally true (noble metals are expensive catalysts), it is not the most relevant or accurate answer for this specific problem. Noble metal catalysts like palladium typically produce hyperbranched structures through a "chain-walking" mechanism, which is different from the formation of "regular branches" via comonomer incorporation. This statement is largely a distractor.

*   **Statement C: "One can use a catalyst of a group VIa transition metal in combination with specific activators."**
    *   **Analysis:** There is an overwhelming consensus among the candidates that this statement is **factually correct and fundamental** to solving the problem. Chromium (Cr), a Group VIa metal, is the key component in the highly selective industrial catalysts used for ethylene trimerization to 1-hexene. This is the cornerstone technology for the "essential additional reaction step."

*   **Statement D: "Aluminum-based activators do not work for the essential additional reaction step."**
    *   **Analysis:** There is an overwhelming consensus that this statement is **factually incorrect**. The chromium-based catalysts described for statement C are almost always activated by aluminum-based compounds (co-catalysts) like methylaluminoxane (MAO), MMAO, or other organoaluminum reagents. These activators are essential for the catalyst to function.

**Step 3: Synthesize and Conclude**
Based on the analysis of the reasoning from all candidates:
*   Statement D is definitively incorrect.
*   Statement B is a plausible but irrelevant distraction.
*   Statement A is ambiguous and its correctness depends on a strict interpretation of "combined system."
*   Statement C is a direct, unambiguous, and fundamentally correct description of the key chemical technology required.

Therefore, the most accurate and reliable statement among the four choices is C. It correctly identifies the class of catalyst that is the industrial and scientific workhorse for the selective oligomerization of ethylene, which is the critical step in forming a polymer with regular branches from a single ethylene feed.

<<<C>>>
"""

# Run the check
result = check_correctness_of_answer(question, llm_answer)
print(result)