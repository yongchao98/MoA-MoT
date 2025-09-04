def check_smeft_symmetry_answer():
    """
    Checks the correctness of an answer about required symmetries in SMEFT.

    The function encodes the fundamental principles of SMEFT construction and
    compares the provided answer against them.
    """

    # 1. Define the established physics principles for SMEFT symmetries.
    # This dictionary represents the ground truth for the check.
    smeft_principles = {
        "Lorentz Symmetry": {
            "is_required": True,
            "justification": "SMEFT is a relativistic quantum field theory, and Lorentz invariance is a foundational requirement for any such theory."
        },
        "Poincare symmetry": {
            "is_required": True,
            "justification": "Poincare symmetry (Lorentz + translations) is the fundamental spacetime symmetry of local QFTs. All operators in the SMEFT Lagrangian must be invariant under it."
        },
        "CP symmetry": {
            "is_required": False,
            "justification": "CP symmetry is known to be violated in the Standard Model's weak sector. Therefore, it is not imposed on SMEFT operators, which can be new sources of CP violation."
        },
        "CPT symmetry": {
            "is_required": True,
            "justification": "The CPT theorem, a robust result in QFT, states that any local, Lorentz-invariant theory must be CPT invariant. SMEFT satisfies these conditions."
        }
    }

    # 2. Map the question's numeric labels to the symmetry names.
    symmetry_map = {
        1: "Lorentz Symmetry",
        2: "Poincare symmetry",
        3: "CP symmetry",
        4: "CPT symmetry"
    }

    # 3. Define the options from the question.
    options = {
        "A": {3, 4},
        "B": {1, 2},
        "C": {1, 3, 4},
        "D": {1, 2, 4}
    }

    # 4. The answer provided by the LLM to be checked.
    llm_answer_choice = "D"

    # 5. Determine the set of required symmetries based on the encoded principles.
    required_symmetries = {
        idx for idx, name in symmetry_map.items()
        if smeft_principles[name]["is_required"]
    }

    # 6. Get the set of symmetries corresponding to the LLM's answer.
    llm_answered_symmetries = options.get(llm_answer_choice)

    if llm_answered_symmetries is None:
        return f"Error: The provided answer '{llm_answer_choice}' is not a valid option."

    # 7. Compare the LLM's answer with the correct set and generate a result.
    if llm_answered_symmetries == required_symmetries:
        return "Correct"
    else:
        # Identify what is wrong with the answer.
        missing = required_symmetries - llm_answered_symmetries
        incorrectly_included = llm_answered_symmetries - required_symmetries
        
        error_messages = []
        if missing:
            for idx in sorted(list(missing)):
                name = symmetry_map[idx]
                justification = smeft_principles[name]["justification"]
                error_messages.append(
                    f"The answer is incorrect because it fails to include the required symmetry: '{name}' ({idx}). "
                    f"Reason: {justification}"
                )
        
        if incorrectly_included:
            for idx in sorted(list(incorrectly_included)):
                name = symmetry_map[idx]
                justification = smeft_principles[name]["justification"]
                error_messages.append(
                    f"The answer is incorrect because it includes a symmetry that is not required: '{name}' ({idx}). "
                    f"Reason: {justification}"
                )
        
        return "\n".join(error_messages)

# Execute the check and print the result.
result = check_smeft_symmetry_answer()
print(result)