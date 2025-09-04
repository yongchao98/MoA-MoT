import re

def check_answer_correctness(question_text, llm_answer_text):
    """
    Checks the correctness of an LLM's answer to the Maxwell's equations question.

    Args:
        question_text (str): The text of the question including the options.
        llm_answer_text (str): The full text of the LLM's answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the error.
    """
    # --- Step 1: Define the ground truth based on physics ---
    # The existence of magnetic monopoles (magnetic charges) means the divergence of B is no longer zero.
    # This modifies Gauss's Law for Magnetism.
    # The existence of magnetic currents (moving monopoles) means the curl of E must be modified for symmetry.
    # This modifies Faraday's Law of Induction.
    correct_changes = {"divergence of the magnetic field", "circulation of the electric field"}

    # --- Step 2: Parse the options from the question text ---
    # This makes the checker robust even if the option letters change.
    options = {}
    # Using a more robust regex to capture the option letter and its text
    option_matches = re.findall(r'([A-D])\)\s(.*?)\n', question_text)
    
    if not option_matches:
        return "Error: Could not parse the multiple-choice options from the question text."

    for letter, text in option_matches:
        # Normalize the text for easier matching
        text = text.lower().strip().replace('.', '')
        claims = set()
        if "divergence of the magnetic field" in text or "divergence and the curl of the magnetic field" in text:
            claims.add("divergence of the magnetic field")
        if "circulation of the electric field" in text:
            claims.add("circulation of the electric field")
        if "circulation of the magnetic field" in text or "curl of the magnetic field" in text:
            claims.add("circulation of the magnetic field")
        if "flux of the electric field" in text:
            claims.add("flux of the electric field")
        options[letter] = claims

    # --- Step 3: Extract the chosen answer from the LLM's response ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The final answer is not in the required format '<<<A>>>', '<<<B>>>', etc."

    chosen_option_letter = match.group(1)
    
    if chosen_option_letter not in options:
        return f"Incorrect. The chosen option '{chosen_option_letter}' is not a valid option (A, B, C, or D)."

    chosen_option_claims = options[chosen_option_letter]

    # --- Step 4: Evaluate the chosen answer against the ground truth ---
    # The chosen option must identify ALL the correct changes and NO incorrect changes.
    # This means the set of claims must be identical to the set of correct changes.
    
    # Check for missing changes
    missing_changes = correct_changes - chosen_option_claims
    if missing_changes:
        return f"Incorrect. The answer is incomplete. It fails to mention that the equation(s) related to '{', '.join(missing_changes)}' would also change."

    # Check for incorrect claims
    incorrect_claims = chosen_option_claims - correct_changes
    if incorrect_claims:
        return f"Incorrect. The answer claims that the equation(s) related to '{', '.join(incorrect_claims)}' would change, which is not true."

    # If both checks pass, the sets are identical.
    return "Correct"

# The question and the final answer provided by the user
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
In a parallel universe where a magnet can have an isolated North or South pole, Maxwell’s equations look different. But, specifically, which of those equations are different?

A) The one related to the divergence of the magnetic field.
B) The one related to the circulation of the magnetic field and the flux of the electric field.
C) The ones related to the divergence and the curl of the magnetic field.
D) The ones related to the circulation of the electric field and the divergence of the magnetic field.

Once you feel you are ready for the final answer, directly return the answer with the format <<<answer content>>> at the end of your response, e.g. <<<C>>>, <<<A>>>
"""

llm_answer = """
Here is a step-by-step analysis of the question, incorporating a review of the candidate answers.

### Step-by-step Analysis of the Physics

To determine the correct answer, we must analyze the physical implications of magnetic monopoles on the four fundamental equations of electromagnetism.

**Careful Point 1: Understanding Standard Maxwell's Equations**

First, let's recall the standard Maxwell's equations and what they describe:
1.  **Gauss's Law for Electricity (∇ ⋅ E = ρ/ε₀):** This relates the **divergence** of the electric field to its source, the electric charge (ρ). It essentially states that electric field lines begin and end on charges.
2.  **Gauss's Law for Magnetism (∇ ⋅ B = 0):** This states that the **divergence** of the magnetic field is zero. This is the precise mathematical statement that there are no magnetic monopoles (no magnetic "charges"). Magnetic field lines always form closed loops.
3.  **Faraday's Law of Induction (∇ × E = -∂B/∂t):** This relates the **circulation (or curl)** of the electric field to a time-varying magnetic field.
4.  **Ampère-Maxwell Law (∇ × B = μ₀J + μ₀ε₀∂E/∂t):** This relates the **circulation (or curl)** of the magnetic field to its sources: electric currents (J) and time-varying electric fields.

**Careful Point 2: The First and Most Direct Change**

The existence of a magnetic monopole means there is a "magnetic charge" that acts as a source or sink for the magnetic field.
*   This directly contradicts Gauss's Law for Magnetism (`∇ ⋅ B = 0`), which is a statement about the *absence* of such sources.
*   Therefore, this law must be modified to include a magnetic charge density term (ρ_m), becoming something like `∇ ⋅ B = μ₀ρ_m`.
*   **Conclusion 1:** The equation related to the **divergence of the magnetic field** must change.

**Careful Point 3: The Second Change due to Symmetry**

If magnetic charges exist, they can also move. A moving magnetic charge would constitute a "magnetic current" (J_m).
*   Maxwell's equations have a strong symmetry. The Ampère-Maxwell law shows that an electric current (J) creates a circulating magnetic field (`∇ × B`).
*   To maintain this symmetry, a magnetic current (J_m) should create a circulating electric field (`∇ × E`).
*   This requires adding a new term to Faraday's Law of Induction, modifying it to something like `∇ × E = -∂B/∂t - μ₀J_m`.
*   **Conclusion 2:** The equation related to the **circulation of the electric field** must also change.

**Careful Point 4: Evaluating the Options**

Based on the analysis, the two equations that must be modified are the one governing the **divergence of the magnetic field** and the one governing the **circulation of the electric field**. Now we match this conclusion to the options provided in the question:

*   A) The one related to the divergence of the magnetic field. — *This is true but incomplete, as a second equation also changes.*
*   B) The one related to the circulation of the magnetic field and the flux of the electric field. — *Incorrect. These are the Ampère-Maxwell Law and Gauss's Law for Electricity, which are not modified.*
*   C) The ones related to the divergence and the curl of the magnetic field. — *Incorrect. The curl of the magnetic field (Ampère-Maxwell Law) is not modified.*
*   D) The ones related to the circulation of the electric field and the divergence of the magnetic field. — **This correctly and completely identifies both required changes.**

### Analysis of the Candidate Answers

A review of the 14 candidate answers reveals a crucial pattern:
*   **Universal Consensus on Reasoning:** Every single candidate answer, without exception, correctly identifies the two physical laws that would change: Gauss's Law for Magnetism (related to the **divergence of the magnetic field**) and Faraday's Law of Induction (related to the **circulation of the electric field**).
*   **Discrepancy in Final Letter:** The final letter chosen (A, B, C, or D) varies significantly between the answers. For example, Answers 2, 6, 7, and 9 all provide the same correct reasoning but select different final letters (A, C, B, and D, respectively). This indicates that the multiple-choice options were presented in a different, randomized order to each agent.

Therefore, the correct approach is to trust the unanimous and correct physical reasoning shared by all candidates and map it to the specific options provided in this prompt. The consensus reasoning points directly to the modification of the laws for the **circulation of the electric field** and the **divergence of the magnetic field**.

<<<D>>>
"""

# Run the check
result = check_answer_correctness(question, llm_answer)
print(result)