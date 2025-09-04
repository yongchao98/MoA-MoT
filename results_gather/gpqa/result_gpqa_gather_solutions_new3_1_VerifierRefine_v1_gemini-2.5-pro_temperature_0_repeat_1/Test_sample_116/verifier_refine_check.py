import re

def check_optical_isomerism_answer(llm_output_text):
    """
    Checks the correctness of the LLM's answer for the optical isomerism question.

    The function verifies if the selected option correctly identifies the compounds
    that exhibit optical isomerism based on established chemical principles.

    Args:
        llm_output_text (str): The full text of the LLM's response, which must
                               contain an answer in the format <<<X>>>.

    Returns:
        str: "Correct" if the answer is right, otherwise a string explaining the error.
    """
    # Step 1: Define the ground truth based on chemical principles.
    # A compound is optically active if it is chiral (non-superimposable on its mirror image).
    compounds_properties = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_chiral": True,
            "reason": "it exhibits atropisomerism (a form of axial chirality) due to restricted rotation around the biphenyl bond caused by bulky ortho substituents. This makes the molecule non-superimposable on its mirror image."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_chiral": True,
            "reason": "it contains a chiral center. The carbon at position 2 is bonded to four different groups (-H, -OH, -CH3, and -COOCH3)."
        },
        3: {
            "name": "benzophenone",
            "is_chiral": False,
            "reason": "it is achiral. It possesses a plane of symmetry that passes through the C=O bond, making it superimposable on its mirror image."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_chiral": False,
            "reason": "it is achiral. The molecule is planar and has a center of inversion, making it superimposable on its mirror image."
        }
    }

    # The set of indices for the correct optically active compounds.
    correct_indices = {idx for idx, prop in compounds_properties.items() if prop["is_chiral"]}
    # Expected correct set: {1, 2}

    # Step 2: Define the mapping from question options to compound indices.
    options_map = {
        "A": {3, 4},
        "B": {1, 2, 4},
        "C": {1, 2},
        "D": {2, 3}
    }

    # Step 3: Parse the LLM's chosen option from its output.
    match = re.search(r'<<<([A-D])>>>', llm_output_text)
    if not match:
        return "Error: Could not find a valid answer in the format <<<X>>> in the provided text."

    llm_choice = match.group(1)
    llm_selected_indices = options_map.get(llm_choice)

    if llm_selected_indices is None:
        return f"Error: The parsed choice '{llm_choice}' is not a valid option (A, B, C, or D)."

    # Step 4: Compare the LLM's answer with the correct answer and generate a report.
    if llm_selected_indices == correct_indices:
        return "Correct"
    else:
        # Build a detailed error message.
        error_message = f"Incorrect. The provided answer chose option {llm_choice}, which corresponds to compounds {sorted(list(llm_selected_indices))}.\n"
        error_message += f"The correct answer is C, which corresponds to compounds {sorted(list(correct_indices))}.\n\n"

        # Explain what was missed.
        missed_compounds = correct_indices - llm_selected_indices
        if missed_compounds:
            error_message += "The answer failed to identify the following optically active compound(s):\n"
            for idx in sorted(list(missed_compounds)):
                error_message += f" - Compound {idx} ({compounds_properties[idx]['name']}) is optically active because {compounds_properties[idx]['reason']}\n"

        # Explain what was incorrectly included.
        incorrectly_included_compounds = llm_selected_indices - correct_indices
        if incorrectly_included_compounds:
            error_message += "\nThe answer incorrectly identified the following compound(s) as optically active:\n"
            for idx in sorted(list(incorrectly_included_compounds)):
                error_message += f" - Compound {idx} ({compounds_properties[idx]['name']}) is NOT optically active because {compounds_properties[idx]['reason']}\n"

        return error_message.strip()

# The final answer provided by the user to be checked.
final_answer_from_user = """
An organic compound shows optical isomerism if it is chiral, meaning its molecule is non-superimposable on its mirror image. This property typically arises from the presence of a chiral center or from restricted rotation that creates a chiral axis (atropisomerism). Let's analyze each compound step-by-step.

**1. dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate**
*   **Structure:** This is a biphenyl molecule, which consists of two phenyl rings connected by a single C-C bond.
*   **Analysis:** The molecule has bulky substituents at all four *ortho* positions (positions 2, 2', 6, and 6'). The substituents are methoxycarbonyl groups (-COOCH₃) and nitro groups (-NO₂). These large groups cause significant steric hindrance, which prevents the two phenyl rings from rotating freely around the single bond that connects them.
*   **Careful Point (Atropisomerism):** This restricted rotation locks the molecule into a twisted, non-planar conformation. Because the substitution on each ring is not symmetric (e.g., one ring has -COOCH₃ at position 2 and -NO₂ at position 6), the molecule as a whole lacks a plane of symmetry or a center of inversion.
*   **Conclusion:** The molecule is chiral due to this axial chirality (a phenomenon called atropisomerism) and is therefore optically active.

**2. methyl 2-hydroxypropanoate**
*   **Structure:** The chemical structure is CH₃-CH(OH)-COOCH₃. It is also known as methyl lactate.
*   **Analysis:** We need to check for a chiral center, which is a carbon atom bonded to four different groups.
*   **Careful Point (Chiral Center):** Let's examine the carbon atom at position 2 (the one attached to the -OH group). The four groups attached to it are:
    1.  A hydrogen atom (-H)
    2.  A hydroxyl group (-OH)
    3.  A methyl group (-CH₃)
    4.  A methoxycarbonyl group (-COOCH₃)
*   **Conclusion:** Since all four groups are different, this carbon is a chiral center. A molecule with a single chiral center is always chiral and thus exhibits optical isomerism.

**3. benzophenone**
*   **Structure:** The structure is (C₆H₅)₂C=O, which is a central carbonyl group attached to two phenyl rings.
*   **Analysis:** The molecule does not have any chiral centers. The central carbon is sp² hybridized.
*   **Careful Point (Symmetry):** The molecule possesses a plane of symmetry that passes through the C=O bond and bisects the angle between the two phenyl rings.
*   **Conclusion:** Due to the presence of a plane of symmetry, the molecule is achiral and superimposable on its mirror image. It does not show optical isomerism.

**4. dimethyl fumarate**
*   **Structure:** This is the dimethyl ester of fumaric acid, which is the *trans*-isomer of but-2-enedioic acid. Its structure is CH₃OOC-CH=CH-COOCH₃, with the two ester groups on opposite sides of the C=C double bond.
*   **Analysis:** The molecule is planar.
*   **Careful Point (Symmetry):** A planar molecule has the plane of the molecule itself as a plane of symmetry. Furthermore, this molecule has a center of inversion (a point of symmetry) at the midpoint of the C=C double bond.
*   **Conclusion:** The presence of these symmetry elements makes the molecule achiral. Therefore, it does not show optical isomerism.

**Final Decision:**
Based on the analysis, only compounds 1 and 2 are chiral and exhibit optical isomerism. Compound 1 shows atropisomerism (axial chirality), and compound 2 has a chiral center. Compounds 3 and 4 are achiral due to the presence of symmetry elements. Therefore, the correct option is the one that includes 1 and 2.

<<<C>>>
"""

# Execute the check
result = check_optical_isomerism_answer(final_answer_from_user)
print(result)