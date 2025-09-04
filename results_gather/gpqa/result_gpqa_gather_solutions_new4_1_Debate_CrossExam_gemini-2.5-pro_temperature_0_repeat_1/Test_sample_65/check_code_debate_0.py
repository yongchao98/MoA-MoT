import re

def check_answer(question_text, llm_answer_text):
    """
    Checks the correctness of an LLM's answer to a multiple-choice question about SMEFT symmetries.

    The function encodes the known physical principles of SMEFT as a "ground truth" and compares
    the selected option against it.

    Args:
        question_text (str): The text of the question, including the options.
        llm_answer_text (str): The full text of the LLM's answer, ending with the choice like <<<B>>>.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the reason for the error.
    """

    # --- Ground Truth Definition ---
    # Based on fundamental principles of Quantum Field Theory and SMEFT construction:
    # 1. Lorentz Symmetry: Required. SMEFT is a relativistic QFT.
    # 2. Poincare Symmetry: Required. Includes Lorentz symmetry and spacetime translations.
    # 3. CP Symmetry: Not required. The SM violates it, and SMEFT parameterizes new sources of violation.
    # 4. CPT Symmetry: Required. A consequence of the CPT theorem for local, Lorentz-invariant QFTs.
    
    symmetries_database = {
        1: {"name": "Lorentz Symmetry", "required": True},
        2: {"name": "Poincare symmetry", "required": True},
        3: {"name": "CP symmetry", "required": False},
        4: {"name": "CPT symmetry", "required": True}
    }
    
    correct_symmetries_set = {num for num, data in symmetries_database.items() if data["required"]}

    # --- Parse Options from Question ---
    # This makes the checker robust even if the options change.
    options = {}
    option_lines = re.findall(r"[A-D]\)\s*([\d,]+)", question_text)
    if not option_lines:
        return "Error: Could not parse the multiple-choice options (A, B, C, D) from the question text."
        
    for i, line in enumerate(option_lines):
        key = chr(ord('A') + i)
        try:
            nums = {int(n.strip()) for n in line.split(',')}
            options[key] = nums
        except ValueError:
            return f"Error: Could not parse numbers for option {key}."

    # --- Parse the LLM's final answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer does not contain a final choice in the format <<<A>>>, <<<B>>>, etc."

    selected_option_key = match.group(1)
    
    if selected_option_key not in options:
        return f"Incorrect: The selected option '{selected_option_key}' is not a valid choice."

    selected_symmetries_set = options[selected_option_key]

    # --- Compare and Generate Result ---
    if selected_symmetries_set == correct_symmetries_set:
        return "Correct"
    else:
        # Generate a detailed reason for the incorrect answer.
        reason = f"The final answer '{selected_option_key}' is incorrect.\n"
        
        missing_symmetries = correct_symmetries_set - selected_symmetries_set
        if missing_symmetries:
            missing_names = [symmetries_database[num]["name"] for num in sorted(list(missing_symmetries))]
            reason += f"Constraint not satisfied: The answer fails to include the required symmetry/symmetries: {', '.join(missing_names)}.\n"
            
        incorrectly_included = selected_symmetries_set - correct_symmetries_set
        if incorrectly_included:
            incorrect_names = [symmetries_database[num]["name"] for num in sorted(list(incorrectly_included))]
            reason += f"Constraint not satisfied: The answer incorrectly includes the symmetry/symmetries which are not required for all operators: {', '.join(incorrect_names)}.\n"
            if 3 in incorrectly_included:
                reason += "Specifically, CP symmetry is not required because the Standard Model itself violates it, and SMEFT is a framework to parameterize new sources of CP violation."

        return reason.strip()

# --- Provided data from the user ---
question = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
Which of the following symmetries must be respected by all operators in the SMEFT?

1. Lorentz Symmetry
2. Poincare symmetry
3. CP symmetry
4. CPT symmetry

A) 3,4
B) 1,2,4
C) 1,3,4
D) 1,2
"""

llm_answer = """
Here is a step-by-step analysis to determine the correct answer.

**Step 1: Understand the Standard Model Effective Field Theory (SMEFT)**
The Standard Model Effective Field Theory (SMEFT) is a theoretical framework that extends the Standard Model (SM). It assumes that any new physics exists at an energy scale (Λ) much higher than the energies accessible by current experiments. The effects of this new physics at lower energies are described by adding higher-dimensional operators (operators with a mass dimension greater than 4) to the Standard Model Lagrangian. A crucial point is that these new operators must be constructed using only SM fields and must respect the fundamental symmetries of the underlying theoretical framework.

**Step 2: Analyze each symmetry individually.**

*   **1. Lorentz Symmetry:** This is a fundamental spacetime symmetry of special relativity, which states that the laws of physics are the same for all observers in uniform motion (i.e., invariant under rotations and boosts). The Standard Model is a relativistic quantum field theory (QFT), and Lorentz invariance is a foundational principle. SMEFT is built upon this same foundation, so all operators in its Lagrangian must be constructed to be Lorentz invariant (specifically, Lorentz scalars). Therefore, Lorentz symmetry **must** be respected.

*   **2. Poincare Symmetry:** This is the full group of spacetime symmetries in special relativity. It includes Lorentz transformations (rotations and boosts) and spacetime translations. Invariance under translations ensures the conservation of energy and momentum. Like Lorentz symmetry, Poincare symmetry is a foundational assumption for any local, relativistic QFT, including the SM and SMEFT. Therefore, Poincare symmetry **must** be respected.

*   **3. CP Symmetry:** This is the symmetry under the combined operation of Charge Conjugation (C) and Parity (P). It is a well-established experimental fact that the weak interactions within the Standard Model violate CP symmetry. Since SMEFT is an extension of the SM, it is not required to be CP-symmetric. In fact, one of the primary goals of the SMEFT program is to parameterize and search for new potential sources of CP violation from physics beyond the Standard Model. Thus, CP symmetry is **not** a required symmetry for all operators.

*   **4. CPT Symmetry:** This is the symmetry under the combined operation of Charge Conjugation (C), Parity (P), and Time Reversal (T). The **CPT Theorem** is a fundamental result in quantum field theory. It states that any local, Lorentz-invariant QFT with a Hermitian Hamiltonian must be CPT invariant. Since SMEFT is constructed to be a local and Lorentz-invariant theory, it must respect CPT symmetry.

**Step 3: Synthesize the results and select the correct option.**

Based on the analysis:
*   Lorentz Symmetry (1) **is required**.
*   Poincare Symmetry (2) **is required**.
*   CP Symmetry (3) is **not required**.
*   CPT Symmetry (4) **is required**.

The set of symmetries that must be respected by all operators in the SMEFT is {1, 2, 4}.

This corresponds to option B.

<<<B>>>
"""

# Execute the check
result = check_answer(question, llm_answer)
print(result)