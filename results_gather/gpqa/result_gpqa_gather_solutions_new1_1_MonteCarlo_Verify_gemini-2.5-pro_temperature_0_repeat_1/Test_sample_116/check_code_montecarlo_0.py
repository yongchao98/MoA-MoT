import re

def check_correctness_of_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The question asks to identify which of the following compounds show optical isomerism:
    1. dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate
    2. methyl 2-hydroxypropanoate
    3. benzophenone
    4. dimethyl fumarate

    Options:
    A) 2 and 3
    B) 1, 2 and 4
    C) 3 and 4
    D) 1 and 2

    Args:
        llm_answer_text: The full text of the LLM's answer, including the final choice in <<<X>>> format.

    Returns:
        "Correct" if the answer is correct, otherwise a string explaining the error.
    """

    # Step 1: Define the ground truth based on chemical principles.
    # A compound shows optical isomerism if it is chiral (non-superimposable on its mirror image).
    compound_analysis = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_optically_active": True,
            "reason": "It exhibits atropisomerism (a form of axial chirality) due to hindered rotation around the central C-C bond caused by bulky ortho substituents. This makes the molecule chiral."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_optically_active": True,
            "reason": "It has a chiral center. The carbon at position 2 is bonded to four different groups (-H, -OH, -CH3, -COOCH3), making the molecule chiral."
        },
        3: {
            "name": "benzophenone",
            "is_optically_active": False,
            "reason": "It is achiral. The molecule has a plane of symmetry that bisects the C=O bond and the two phenyl rings."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_optically_active": False,
            "reason": "It is achiral. The molecule is planar and has a center of inversion, making it superimposable on its mirror image."
        }
    }

    # Determine the set of compounds that are actually optically active.
    correct_set = {
        compound_id for compound_id, properties in compound_analysis.items()
        if properties["is_optically_active"]
    }

    # Step 2: Define the sets corresponding to each multiple-choice option.
    options = {
        "A": {2, 3},
        "B": {1, 2, 4},
        "C": {3, 4},
        "D": {1, 2}
    }

    # Step 3: Parse the final answer from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The final answer is not provided in the required format '<<<X>>>' where X is one of A, B, C, or D."

    llm_choice = match.group(1)
    llm_set = options.get(llm_choice)

    if llm_set is None:
        # This case should not happen with the regex used, but it's good practice.
        return f"Incorrect. The answer choice '{llm_choice}' is not a valid option."

    # Step 4: Compare the LLM's answer with the ground truth.
    if llm_set == correct_set:
        return "Correct"
    else:
        # Find the correct option letter for the error message.
        correct_option_letter = [key for key, val in options.items() if val == correct_set][0]
        
        error_message = (
            f"Incorrect. The provided answer is '{llm_choice}', which corresponds to compounds {sorted(list(llm_set))}.\n"
            f"The correct answer is '{correct_option_letter}', which corresponds to compounds {sorted(list(correct_set))}.\n\n"
            "Reasoning:\n"
            "A compound shows optical isomerism if it is chiral (lacks a plane of symmetry or a center of inversion).\n"
        )

        # Provide a breakdown of why each compound is or isn't optically active.
        for i in range(1, 5):
            is_active = compound_analysis[i]["is_optically_active"]
            reason = compound_analysis[i]["reason"]
            
            error_message += f"- Compound {i}: Is {'optically active' if is_active else 'NOT optically active'}. {reason}\n"

        # Highlight the specific mistake in the LLM's choice.
        missed = correct_set - llm_set
        extra = llm_set - correct_set
        
        if missed:
            error_message += f"\nThe answer '{llm_choice}' incorrectly omits optically active compound(s): {sorted(list(missed))}.\n"
        if extra:
            error_message += f"The answer '{llm_choice}' incorrectly includes achiral compound(s): {sorted(list(extra))}.\n"
            
        return error_message

# The final answer provided by the user to be checked.
llm_answer_to_check = """
To determine which compounds show optical isomerism, we must assess each molecule for chirality. A molecule is chiral if it is non-superimposable on its mirror image. The presence of a chiral center or an element of axial chirality (like in atropisomers) are common reasons for a molecule to be chiral, provided it does not have an internal plane of symmetry or a center of inversion.

Here is a step-by-step analysis of each compound:

1.  **dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate**:
    *   **Structure**: This is a biphenyl system, where two benzene rings are connected by a single C-C bond.
    *   **Careful Point (Atropisomerism)**: The key feature is the presence of bulky substituents (nitro groups and dimethyl carboxylate groups) at all four *ortho* positions (2, 2', 6, and 6'). These large groups sterically hinder free rotation around the single bond connecting the rings.
    *   **Analysis**: This restricted rotation locks the molecule into a stable, twisted, non-planar conformation. This conformation lacks both a plane of symmetry and a center of inversion, making the molecule chiral. This specific type of chirality arising from hindered rotation is called atropisomerism.
    *   **Conclusion**: The compound is chiral and therefore exhibits optical isomerism.

2.  **methyl 2-hydroxypropanoate**:
    *   **Structure**: Also known as methyl lactate, its formula is CH₃-CH(OH)-COOCH₃.
    *   **Careful Point (Chiral Center)**: We must check for carbon atoms bonded to four different groups.
    *   **Analysis**: The carbon atom at position 2 is bonded to four distinct groups:
        1.  -H (a hydrogen atom)
        2.  -OH (a hydroxyl group)
        3.  -CH₃ (a methyl group)
        4.  -COOCH₃ (a methoxycarbonyl group)
    *   **Conclusion**: Since the carbon at position 2 is a chiral center, the molecule is chiral and exhibits optical isomerism.

3.  **benzophenone**:
    *   **Structure**: (C₆H₅)₂C=O.
    *   **Careful Point (Symmetry)**: We must check for elements of symmetry.
    *   **Analysis**: The molecule has no chiral centers. It possesses a plane of symmetry that passes through the C=O bond and bisects the angle between the two phenyl rings. Any molecule with a plane of symmetry is achiral.
    *   **Conclusion**: The compound is achiral and does not show optical isomerism.

4.  **dimethyl fumarate**:
    *   **Structure**: This is the *trans*-isomer of dimethyl but-2-enedioate.
    *   **Careful Point (Symmetry)**: This is a planar molecule.
    *   **Analysis**: The molecule is planar and possesses a center of inversion at the midpoint of the C=C double bond. The presence of a center of inversion (or a plane of symmetry) makes a molecule achiral.
    *   **Conclusion**: The compound is achiral and does not show optical isomerism.

**Final Decision**:
Based on the analysis, compounds 1 and 2 are chiral and will show optical isomerism. Compounds 3 and 4 are achiral and will not. Therefore, the correct option includes only 1 and 2.

This corresponds to option D.

<<<D>>>
"""

# Execute the check and print the result.
result = check_correctness_of_answer(llm_answer_to_check)
print(result)