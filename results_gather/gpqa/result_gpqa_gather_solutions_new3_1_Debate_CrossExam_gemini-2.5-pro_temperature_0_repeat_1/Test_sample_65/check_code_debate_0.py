import re

def check_correctness(question, llm_answer):
    """
    Checks the correctness of the LLM's answer for the SMEFT symmetries question.

    Args:
        question (str): The question text.
        llm_answer (str): The full response from the LLM, including the final answer in <<<>>> format.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # --- Ground Truth Definition ---
    # Based on the principles of Quantum Field Theory and the construction of SMEFT:
    # 1. Lorentz Symmetry: Required. Foundational for relativistic QFT.
    # 2. Poincare Symmetry: Required. The full spacetime symmetry group for special relativity.
    # 3. CP Symmetry: Not required. The SM itself violates CP, and SMEFT parameterizes new sources of CP violation.
    # 4. CPT Symmetry: Required. A consequence of the CPT theorem for local, Lorentz-invariant QFTs.
    
    correct_symmetries = {1, 2, 4}
    symmetries_map = {
        1: "Lorentz Symmetry",
        2: "Poincare symmetry",
        3: "CP symmetry",
        4: "CPT symmetry"
    }

    # --- Options Mapping ---
    # A) 3,4
    # B) 1,2,4
    # C) 1,3,4
    # D) 1,2
    options = {
        'A': {3, 4},
        'B': {1, 2, 4},
        'C': {1, 3, 4},
        'D': {1, 2}
    }

    # --- Extract LLM's Final Answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Invalid answer format. The final answer should be in the format <<<A>>>, <<<B>>>, etc."

    chosen_option = match.group(1)
    
    if chosen_option not in options:
        return f"Invalid option '{chosen_option}'. The option must be one of {list(options.keys())}."

    chosen_symmetries = options[chosen_option]

    # --- Verification ---
    if chosen_symmetries == correct_symmetries:
        return "Correct"
    else:
        # Generate a detailed reason for the error
        missing_symmetries = correct_symmetries - chosen_symmetries
        incorrectly_included_symmetries = chosen_symmetries - correct_symmetries
        
        error_messages = []
        if missing_symmetries:
            missing_names = [f"'{symmetries_map[s]} ({s})'" for s in sorted(list(missing_symmetries))]
            error_messages.append(f"it incorrectly omits the required symmetr{'y' if len(missing_names) == 1 else 'ies'}: {', '.join(missing_names)}")
        
        if incorrectly_included_symmetries:
            included_names = [f"'{symmetries_map[s]} ({s})'" for s in sorted(list(incorrectly_included_symmetries))]
            error_messages.append(f"it incorrectly includes the non-required symmetr{'y' if len(included_names) == 1 else 'ies'}: {', '.join(included_names)}")

        reason = " and ".join(error_messages)
        
        return (f"Incorrect. The chosen option '{chosen_option}' is wrong because {reason}. "
                f"The correct set of required symmetries is {{1, 2, 4}} (Lorentz, Poincare, and CPT).")


# The user's provided question and the final answer from the LLM analysis
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
Here is a step-by-step analysis to determine the final answer.

### Step 1: Analyze the Definition of SMEFT
The Standard Model Effective Field Theory (SMEFT) is an extension of the Standard Model (SM). It is constructed by adding higher-dimensional operators to the SM Lagrangian. A core principle of this construction is that these new operators must respect the fundamental symmetries of the underlying theoretical framework, which is a relativistic Quantum Field Theory (QFT).

### Step 2: Evaluate Each Symmetry
1.  **Lorentz Symmetry (1):** The Standard Model is a relativistic QFT, meaning it is built on the principles of special relativity. A cornerstone of special relativity is Lorentz invariance, which dictates that the laws of physics are the same for all inertial observers. SMEFT inherits this foundational requirement. All operators are constructed to be Lorentz scalars, meaning they are invariant under Lorentz transformations. Therefore, Lorentz symmetry **must** be respected.

2.  **Poincare Symmetry (2):** The Poincare group represents the full spacetime symmetries of special relativity. It includes Lorentz transformations (rotations and boosts) and spacetime translations. Invariance under translations ensures the conservation of energy and momentum. The SMEFT Lagrangian is required to be invariant under the Poincare group. Since the Lorentz group is a subgroup of the Poincare group, respecting Poincare symmetry automatically implies respecting Lorentz symmetry. Therefore, Poincare symmetry **must** be respected.

3.  **CP Symmetry (3):** CP symmetry is the combination of Charge Conjugation (C) and Parity (P). It is a well-established experimental fact that the Standard Model itself violates CP symmetry (e.g., in the weak interactions, described by the CKM matrix). Since the base theory (the SM) is not CP-symmetric, there is no reason to impose this symmetry on the new operators in SMEFT. In fact, a primary motivation for SMEFT is to provide a framework for parameterizing potential *new* sources of CP violation from high-energy physics. Therefore, CP symmetry is **not** a required symmetry for all operators.

4.  **CPT Symmetry (4):** CPT symmetry is the combination of Charge Conjugation (C), Parity (P), and Time Reversal (T). The CPT theorem is a fundamental result in QFT, stating that any local, Lorentz-invariant quantum field theory with a Hermitian Hamiltonian must be CPT symmetric. Since SMEFT is constructed as a local and Lorentz-invariant QFT, it must obey the CPT theorem. Therefore, CPT symmetry **must** be respected.

### Step 3: Synthesize the Results and Select the Correct Option
Based on the analysis, the symmetries that must be respected by all operators in SMEFT are:
*   Lorentz Symmetry (1)
*   Poincare Symmetry (2)
*   CPT Symmetry (4)

The symmetry that is not required is CP Symmetry (3).

The correct combination of required symmetries is {1, 2, 4}. Now, we match this set to the given options:
*   A) 3,4
*   B) 1,2,4
*   C) 1,3,4
*   D) 1,2

The set {1, 2, 4} corresponds to option B.

### Step 4: Review of Candidate Answers
The majority of the provided answers (1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15) correctly reason that Lorentz, Poincare, and CPT are required symmetries, while CP is not. Several of these answers then make a mistake in selecting the final letter corresponding to the set {1, 2, 4}. However, the physical reasoning is consistent and points to {1, 2, 4}. Answers 2, 4, 5, 9, 13, and 15 correctly identify B as the final answer.

<<<B>>>
"""

# Run the check
result = check_correctness(question, llm_answer)
print(result)