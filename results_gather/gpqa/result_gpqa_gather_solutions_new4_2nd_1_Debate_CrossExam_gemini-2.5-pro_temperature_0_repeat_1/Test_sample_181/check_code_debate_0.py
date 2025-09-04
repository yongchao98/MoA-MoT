import re

def check_correctness_of_answer(question_and_options, llm_answer_block):
    """
    Checks the correctness of the LLM's answer for the Mott-Gurney equation question.

    The function verifies the chosen option against the fundamental assumptions
    of the Mott-Gurney law for Space-Charge-Limited Current (SCLC).
    """

    # --- Step 1: Define the ground truth based on physics principles ---
    # The Mott-Gurney equation is valid under these conditions:
    # 1. Single-carrier device (unipolar injection)
    # 2. Trap-free material
    # 3. Ohmic contact (no carrier injection barrier)
    # 4. Drift-dominated current (diffusion current is negligible)

    # We can represent these as checks. A correct statement must pass all.
    def evaluate_statement(statement):
        """
        Evaluates a single option statement against the physical principles.
        Returns a tuple (is_correct, reason_for_failure).
        """
        statement = statement.lower() # Case-insensitive matching

        # Check for explicitly incorrect conditions
        if "two-carrier" in statement:
            return (False, "The Mott-Gurney law is for a single-carrier device, not a two-carrier device.")
        if "schottky contact" in statement:
            return (False, "A Schottky contact has an injection barrier, which is contrary to the SCLC requirement of an Ohmic contact (no injection barrier).")
        if "negligible drift current" in statement:
            return (False, "SCLC is a drift-dominated current; it is the diffusion current that is negligible, not the drift current.")

        # Check if all required conditions are present in the statement
        missing_conditions = []
        if "single-carrier" not in statement:
            missing_conditions.append("'single-carrier device'")
        if "trap-free" not in statement:
            # Note: An incorrect option might omit a correct detail.
            # We check this to ensure the chosen *correct* option is complete.
            # Option D is wrong for another reason but also omits this.
            pass # Let's not fail for omission if there isn't an explicit error.
        if "no carrier injection barrier" not in statement and "ohmic contact" not in statement:
            missing_conditions.append("'no carrier injection barrier' or 'Ohmic contact'")
        if "negligible diffusion current" not in statement:
            missing_conditions.append("'negligible diffusion current'")

        # The correct answer must contain all key elements and no incorrect ones.
        # Let's check for the combination that makes a statement correct.
        is_correct = (
            "single-carrier" in statement and
            "trap-free" in statement and
            ("no carrier injection barrier" in statement or "ohmic contact" in statement) and
            "negligible diffusion current" in statement
        )

        if is_correct:
            return (True, "")
        elif not missing_conditions:
            # This case should not be reached if the logic above is sound, but as a fallback:
            return (False, "The statement contains a subtle error or is incomplete.")
        else:
            # This handles cases where a statement is incomplete but doesn't have an explicit error.
            # For this specific question, all incorrect options have explicit errors.
            return (False, f"The statement is incomplete and is missing the condition(s): {', '.join(missing_conditions)}.")


    # --- Step 2: Parse the LLM's answer and the options from the prompt ---
    try:
        # Extract the final answer letter (e.g., 'A') from <<<A>>>
        chosen_option_letter_match = re.search(r'<<<([A-D])>>>', llm_answer_block)
        if not chosen_option_letter_match:
            return "Failure: Could not parse the final answer from the provided text (e.g., '<<<A>>>')."
        chosen_option_letter = chosen_option_letter_match.group(1)

        # Parse the options from the question text
        options_raw = re.findall(r'\n([A-D])\)\s(.*?)(?=\n[A-D]\)|$)', question_and_options, re.DOTALL)
        options = {opt[0]: opt[1].strip() for opt in options_raw}

        if not options or len(options) != 4:
             return "Failure: Could not parse the four options (A, B, C, D) from the question text."

        chosen_statement = options.get(chosen_option_letter)
        if chosen_statement is None:
            return f"Failure: The answer '{chosen_option_letter}' does not correspond to any of the parsed options."

    except Exception as e:
        return f"An error occurred during parsing: {e}"

    # --- Step 3: Evaluate the chosen option and return the result ---
    is_correct, reason = evaluate_statement(chosen_statement)

    if is_correct:
        return "Correct"
    else:
        return f"Incorrect. The chosen answer is option {chosen_option_letter}, which states: '{chosen_statement}'. This is wrong because: {reason}"


# The user's provided context
question_text = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
The Mott-Gurney equation describes the dark current (J) versus voltage (V) behavior of a semiconductor device in the space-charge-limited current (SCLC) regime. The equation can be expressed as
$ J = \frac{9}{8} \epsilon \mu \frac{V^2}{L^3}$
where $\epsilon$ is the dielectric constant, $\mu$ is the charge carrier mobility, and L is the length of the device. Which of the following statements is true about the validity of this equation?

A) The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.
B) The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.
C) The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.
D) The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.
"""

llm_response = """
Here is a step-by-step analysis of the question, synthesizing the reasoning from the provided candidate answers and highlighting the key points for a definitive conclusion.

### Step-by-Step Analysis

1.  **Deconstruct the Core Concept:** The question asks for the conditions under which the Mott-Gurney equation for Space-Charge-Limited Current (SCLC) is valid. The equation, $ J = \frac{9}{8} \epsilon \mu \frac{V^2}{L^3}$, describes a specific physical regime. The validity of this equation depends entirely on the assumptions made during its derivation.

2.  **Identify the Fundamental Assumptions:** A thorough review of the candidate answers and the principles of semiconductor physics reveals a strong consensus on four key assumptions:
    *   **Single-Carrier Device:** The model is unipolar, meaning it assumes only one type of charge carrier (e.g., only electrons or only holes) is mobile and contributes to the current. This is why there is a single mobility term ($\mu$) in the equation.
    *   **Trap-Free Material:** The ideal Mott-Gurney law applies to a perfect insulator or semiconductor with no localized energy states (traps) in the bandgap. If traps were present, they would capture some of the injected carriers, altering the current-voltage relationship.
    *   **Ohmic Contact (No Injection Barrier):** For the current to be limited by the "space charge" in the bulk of the material, the injecting contact must be able to supply an essentially unlimited number of carriers. This requires a contact with no energy barrier to injection. This is the definition of an ideal Ohmic contact in this context. A Schottky contact, which has a barrier, would limit the current by the injection process itself, not by the bulk space charge.
    *   **Drift-Dominated Current (Negligible Diffusion):** The model assumes that the primary transport mechanism is the movement of carriers under the influence of the applied electric field (drift). The current component arising from the concentration gradient of carriers (diffusion) is assumed to be negligible. Therefore, stating that *drift* current is negligible is a direct contradiction of the model's basis.

3.  **Evaluate Each Option Against the Assumptions:** Now, we can systematically evaluate the given options from the prompt against these four established conditions.

    *   **A) The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.**
        *   `trap-free`: **Correct.** (Matches Assumption 2)
        *   `single-carrier device`: **Correct.** (Matches Assumption 1)
        *   `no carrier injection barrier`: **Correct.** (Matches Assumption 3)
        *   `negligible diffusion current`: **Correct.** (Matches Assumption 4)
        *   **Conclusion:** This statement correctly and completely combines all the necessary conditions.

    *   **B) The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.**
        *   `two-carrier device`: **Incorrect.** This contradicts the single-carrier assumption.
        *   **Conclusion:** This statement is incorrect.

    *   **C) The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.**
        *   `negligible drift current`: **Incorrect.** This contradicts the drift-dominated current assumption. The SCLC is fundamentally a drift current.
        *   **Conclusion:** This statement is incorrect.

    *   **D) The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.**
        *   `Schottky contact`: **Incorrect.** This contradicts the no-injection-barrier assumption. A Schottky contact is a barrier contact that limits injection.
        *   **Conclusion:** This statement is incorrect.

### Careful Points and Final Decision

*   The analysis shows that the terms "Ohmic contact" and "no carrier injection barrier" are used synonymously in this context to describe the ideal injecting contact required for SCLC.
*   It is crucial to distinguish between drift and diffusion current. The Mott-Gurney model assumes **drift is dominant** and **diffusion is negligible**.
*   The reasoning across the majority of the provided candidate answers is sound and converges on the same set of physical principles, even though some of them selected a different final letter, likely due to a shuffling of the options in their specific prompts.
*   Based on a rigorous evaluation, only one option withstands scrutiny.

<<<A>>>
"""

# Run the check
result = check_correctness_of_answer(question_text, llm_response)
print(result)