import re

def check_maxwell_equations_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the answer to the Maxwell's equations question.

    The question asks which of Maxwell's equations would change if magnetic monopoles existed.
    The correct answer is that two equations would change:
    1. Gauss's Law for Magnetism (related to the divergence of the magnetic field).
    2. Faraday's Law of Induction (related to the circulation of the electric field).
    """

    # 1. Define the correct physical principles that must change.
    # Note: "curl" is synonymous with "circulation", and "flux" through a closed surface
    # is related to "divergence" by the Divergence Theorem.
    correct_changing_concepts = {
        "circulation of the electric field",
        "divergence of the magnetic field"
    }

    # 2. Define the concepts associated with each option from the question.
    options_map = {
        "A": {"circulation of the magnetic field", "flux of the electric field"},
        "B": {"circulation of the electric field", "divergence of the magnetic field"},
        "C": {"divergence of the magnetic field", "curl of the magnetic field"},
        "D": {"divergence of the magnetic field"}
    }
    
    # Normalize terms for comparison (e.g., curl -> circulation, flux -> divergence)
    normalized_options = {}
    for key, concepts in options_map.items():
        normalized_concepts = set()
        for concept in concepts:
            if "curl" in concept:
                normalized_concepts.add(concept.replace("curl", "circulation"))
            elif "flux" in concept:
                normalized_concepts.add(concept.replace("flux", "divergence"))
            else:
                normalized_concepts.add(concept)
        normalized_options[key] = normalized_concepts


    # 3. Find the letter of the correct option.
    correct_option_letter = None
    for letter, concepts in normalized_options.items():
        # The correct answer must be complete, matching both changing concepts.
        if concepts == correct_changing_concepts:
            correct_option_letter = letter
            break

    if not correct_option_letter:
        # This is a fallback for an issue with the checker's own logic.
        return "Error in checker: Could not identify the correct option among the choices."

    # 4. Extract the LLM's chosen answer from the text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "The answer is incorrect because it does not provide a final answer in the required format, e.g., <<<B>>>."
    
    chosen_option = match.group(1)

    # 5. Compare the LLM's choice with the correct answer.
    if chosen_option == correct_option_letter:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness.
        reason = f"The chosen answer is {chosen_option}, but the correct answer is {correct_option_letter}.\n"
        reason += "The existence of magnetic monopoles would change two of Maxwell's equations:\n"
        reason += "1. Gauss's Law for Magnetism (∇ ⋅ B = 0), which describes the 'divergence of the magnetic field'. It must change to account for magnetic charges.\n"
        reason += "2. Faraday's Law of Induction (∇ × E = -∂B/∂t), which describes the 'circulation of the electric field'. It must change to account for magnetic currents for symmetry.\n"
        
        chosen_concepts = normalized_options.get(chosen_option, set())
        
        if chosen_concepts == {"divergence of the magnetic field"}:
            reason += f"The chosen option {chosen_option} is incorrect because it is incomplete. While the law for the 'divergence of the magnetic field' does change, the law for the 'circulation of the electric field' also changes, and the question asks for the 'equations' (plural)."
        else:
            incorrect_concepts = chosen_concepts - correct_changing_concepts
            if incorrect_concepts:
                 reason += f"The chosen option {chosen_option} is incorrect because it mentions concepts that do not change, such as: '{', '.join(incorrect_concepts)}'."

        return reason

# The user-provided text containing the question and the LLM's answer.
llm_response = """
Here is a step-by-step analysis of the question, incorporating careful points to arrive at the correct answer.

### Step-by-Step Analysis

1.  **Understand the Premise:** The question asks us to consider a hypothetical universe where magnetic monopoles (isolated North or South poles) exist. We need to determine which of Maxwell's four fundamental equations of electromagnetism would be different in this scenario.

2.  **Recall the Standard Maxwell's Equations:** To know what changes, we must first understand the standard equations and their physical meaning.
    *   **Gauss's Law for Electricity (`∇ ⋅ E = ρ/ε₀`):** This relates the **divergence** of the electric field to electric charge density. It means electric field lines originate from and terminate on electric charges. This is sometimes referred to as the law for the flux of the electric field.
    *   **Gauss's Law for Magnetism (`∇ ⋅ B = 0`):** This states that the **divergence** of the magnetic field is zero. This is the precise mathematical statement that there are no magnetic monopoles. It means magnetic field lines always form closed loops.
    *   **Faraday's Law of Induction (`∇ × E = -∂B/∂t`):** This relates the **curl** (or **circulation**) of the electric field to a time-varying magnetic field.
    *   **Ampère-Maxwell Law (`∇ × B = μ₀(J + ε₀∂E/∂t)`):** This relates the **curl** (or **circulation**) of the magnetic field to electric currents and changing electric fields.

3.  **Analyze the Impact of Magnetic Monopoles:**
    *   **Careful Point 1: Magnetic Charge.** A magnetic monopole is effectively a "magnetic charge." Just as an electric charge is a source or sink for electric field lines, a magnetic charge would be a source or sink for magnetic field lines.
    *   **Impact on Gauss's Law for Magnetism:** The equation `∇ ⋅ B = 0` is a direct statement of the *absence* of magnetic monopoles. If they existed, the divergence of the magnetic field would no longer be zero; it would be proportional to the magnetic charge density. Therefore, the equation related to the **divergence of the magnetic field** must change.

    *   **Careful Point 2: Symmetry and Magnetic Current.** The standard equations have a beautiful symmetry, broken only by the lack of magnetic monopoles. The Ampère-Maxwell law shows that a moving electric charge (an electric current) creates a circulating magnetic field. If magnetic monopoles exist, a moving magnetic monopole would constitute a "magnetic current."
    *   **Impact on Faraday's Law of Induction:** To restore the symmetry of the laws, this magnetic current must create a circulating electric field. This effect is not accounted for in the standard Faraday's Law. Therefore, a new term for magnetic current must be added to the equation. This means the equation related to the **circulation of the electric field** must also change.

4.  **Evaluate the Options:**
    Based on the analysis, the two equations that would be different are the ones related to:
    1.  The **divergence of the magnetic field** (Gauss's Law for Magnetism).
    2.  The **circulation of the electric field** (Faraday's Law of Induction).

    Let's match this conclusion with the given options:
    *   A) The one related to the circulation of the magnetic field and the flux of the electric field. — *Incorrect.* These are the Ampère-Maxwell Law and Gauss's Law for Electricity, which do not change.
    *   B) The ones related to the circulation of the electric field and the divergence of the magnetic field. — **Correct.** This option perfectly and completely matches our analysis.
    *   C) The ones related to the divergence and the curl of the magnetic field. — *Incorrect.* The curl of the magnetic field (Ampère-Maxwell Law) does not change.
    *   D) The one related to the divergence of the magnetic field. — *Incomplete.* While this is true, Faraday's Law also changes, and the question asks which "equations" (plural) are different. Option B is the most complete answer.

<<<B>>>
"""

# Run the check
result = check_maxwell_equations_answer(llm_response)
print(result)