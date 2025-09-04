import re

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to the chemistry question.
    It does this by:
    1. Defining the ground truth for the optical activity of each compound.
    2. Defining the multiple-choice options.
    3. Determining the correct option based on the ground truth.
    4. Comparing the correct option with the LLM's provided answer.
    """

    # Ground truth based on chemical principles.
    # A compound is optically active if it is chiral.
    compound_properties = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_optically_active": True,
            "reason": "Exhibits atropisomerism due to restricted rotation, making it chiral."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_optically_active": True,
            "reason": "Contains a chiral center bonded to four different groups."
        },
        3: {
            "name": "benzophenone",
            "is_optically_active": False,
            "reason": "Is achiral due to a plane of symmetry."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_optically_active": False,
            "reason": "Is achiral due to a center of inversion and a plane of symmetry."
        }
    }

    # Define the options as provided in the question.
    options = {
        "A": {3, 4},
        "B": {2, 3},
        "C": {1, 2, 4},
        "D": {1, 2}
    }

    # The final answer provided by the LLM.
    llm_answer_text = """
    Here is a step-by-step analysis to determine the final answer, including all the careful points for each compound.

    **Core Concept: Optical Isomerism and Chirality**
    The fundamental requirement for a compound to show optical isomerism is that it must be **chiral**. A molecule is chiral if it is non-superimposable on its mirror image. To determine this, we look for elements of chirality (like chiral centers or a chiral axis) and check for the absence of key symmetry elements (like a plane of symmetry or a center of inversion).

    **Step-by-step Analysis of Each Compound**

    1.  **dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate**
        *   **Structure**: This molecule consists of two phenyl rings linked by a single C-C bond.
        *   **Careful Point (Atropisomerism)**: The key feature is the presence of large, bulky substituents (-NO₂ and -COOCH₃) at all four *ortho*-positions (the positions immediately adjacent to the bond connecting the rings: 2, 2', 6, and 6').
        *   **Analysis**: These bulky groups cause significant steric hindrance, which prevents the two phenyl rings from rotating freely around the central single bond. This restricted rotation locks the molecule into a stable, twisted, non-planar conformation. This specific type of axial chirality is known as **atropisomerism**.
        *   **Conclusion**: The resulting twisted structure lacks both a plane of symmetry and a center of inversion. Therefore, the molecule is chiral and will exhibit optical isomerism.

    2.  **methyl 2-hydroxypropanoate**
        *   **Structure**: Also known as methyl lactate, its formula is CH₃-CH(OH)-COOCH₃.
        *   **Careful Point (Chiral Center)**: We must check for a carbon atom that is bonded to four different groups.
        *   **Analysis**: The carbon atom at position 2 is bonded to four distinct groups:
            1.  A hydrogen atom (-H)
            2.  A hydroxyl group (-OH)
            3.  A methyl group (-CH₃)
            4.  A methoxycarbonyl group (-COOCH₃)
        *   **Conclusion**: Since this carbon is bonded to four different substituents, it is a **chiral center**. A molecule containing a single chiral center is inherently chiral and therefore optically active.

    3.  **benzophenone**
        *   **Structure**: The structure is (C₆H₅)₂C=O, with a central carbonyl group bonded to two phenyl rings.
        *   **Careful Point (Symmetry)**: We must check for elements of symmetry that would make the molecule achiral.
        *   **Analysis**: The molecule does not have a chiral center. Crucially, it possesses a **plane of symmetry** that passes through the C=O bond and bisects the angle between the two identical phenyl rings.
        *   **Conclusion**: Any molecule with a plane of symmetry is achiral and cannot be optically active.

    4.  **dimethyl fumarate**
        *   **Structure**: This is the *trans*-isomer of an alkene, with the structure CH₃OOC-CH=CH-COOCH₃.
        *   **Careful Point (Symmetry)**: The geometry of the double bond is critical. As the *trans*-isomer, the two ester groups are on opposite sides of the double bond.
        *   **Analysis**: The molecule is planar. It possesses both a **plane of symmetry** (the plane of the molecule itself) and, more definitively, a **center of inversion** located at the midpoint of the C=C double bond.
        *   **Conclusion**: The presence of these symmetry elements makes the molecule achiral and thus not optically active.

    **Final Decision**
    Based on the analysis:
    *   Compound **1** is optically active (due to atropisomerism).
    *   Compound **2** is optically active (due to a chiral center).
    *   Compound **3** is not optically active (achiral).
    *   Compound **4** is not optically active (achiral).

    The compounds that show optical isomerism are 1 and 2. This corresponds to option D.

    <<<D>>>
    """

    # Extract the LLM's chosen option letter
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not parse the LLM's final answer from the text."
    llm_choice = match.group(1)

    # Determine the correct set of optically active compounds
    correct_active_set = {num for num, props in compound_properties.items() if props["is_optically_active"]}

    # Find the correct option letter
    correct_option_letter = None
    for letter, compound_set in options.items():
        if compound_set == correct_active_set:
            correct_option_letter = letter
            break

    # Check if the LLM's reasoning within the text matches the ground truth
    llm_reasoning_set = set()
    if "1 is optically active" in llm_answer_text.replace("**", ""):
        llm_reasoning_set.add(1)
    if "2 is optically active" in llm_answer_text.replace("**", ""):
        llm_reasoning_set.add(2)
    
    # Final validation
    if llm_reasoning_set != correct_active_set:
        return (f"Incorrect. The reasoning in the text is flawed. "
                f"The text identifies {sorted(list(llm_reasoning_set))} as optically active, "
                f"but the correct set is {sorted(list(correct_active_set))}.")

    if llm_choice == correct_option_letter:
        return "Correct"
    else:
        return (f"Incorrect. The final answer choice is wrong. "
                f"The reasoning correctly identifies that compounds {sorted(list(correct_active_set))} are optically active, "
                f"which corresponds to option {correct_option_letter}. "
                f"However, the final answer selected was {llm_choice}.")

# Run the check
result = check_correctness()
print(result)