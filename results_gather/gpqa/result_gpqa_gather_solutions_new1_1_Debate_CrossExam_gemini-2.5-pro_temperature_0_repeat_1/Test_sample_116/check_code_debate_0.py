import re

def check_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of an LLM's answer to a chemistry question about optical isomerism.

    The function verifies the answer by:
    1. Defining the correct chirality for each compound based on established chemical principles.
    2. Determining the correct option based on this analysis.
    3. Extracting the LLM's chosen option from the provided text.
    4. Comparing the LLM's choice with the correct option.

    Args:
        llm_answer_text: A string containing the LLM's full response, including the final answer
                         in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or providing a detailed reason for the incorrectness.
    """
    # Step 1: Define the chemical properties of each compound based on established principles.
    # A compound shows optical isomerism if it is chiral.
    compound_analysis = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_chiral": True,
            "reason": "Exhibits atropisomerism (a form of axial chirality) due to bulky ortho substituents restricting rotation around the biphenyl bond. This makes the molecule non-superimposable on its mirror image."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_chiral": True,
            "reason": "Contains a chiral center (the second carbon) which is bonded to four different groups (-H, -OH, -CH3, and -COOCH3)."
        },
        3: {
            "name": "benzophenone",
            "is_chiral": False,
            "reason": "Is achiral because it possesses a plane of symmetry that passes through the C=O bond and bisects the angle between the two phenyl rings."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_chiral": False,
            "reason": "Is achiral because it is a planar molecule with a center of inversion at the midpoint of the C=C double bond."
        }
    }

    # Step 2: Identify the correct set of optically active compounds.
    correct_compounds = {num for num, props in compound_analysis.items() if props["is_chiral"]}
    
    # Step 3: Define the options from the question.
    options = {
        'A': {1, 2},
        'B': {1, 2, 4},
        'C': {2, 3},
        'D': {3, 4}
    }

    # Step 4: Determine the correct option letter.
    correct_option = None
    for option, compounds in options.items():
        if compounds == correct_compounds:
            correct_option = option
            break
            
    if correct_option is None:
        # This case should not be reached if the question is well-formed.
        return "Error in checker: Could not find a matching option for the correct answer."

    # Step 5: Extract the final answer from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The final answer is not provided in the required format '<<<A>>>', '<<<B>>>', etc."
        
    llm_option = match.group(1)

    # Step 6: Compare the LLM's answer with the correct answer and provide feedback.
    if llm_option == correct_option:
        return "Correct"
    else:
        reason_str = f"Incorrect: The provided answer '{llm_option}' is wrong. The correct answer is '{correct_option}'.\n\n"
        reason_str += "Here is the detailed analysis of each compound:\n"
        for num, props in compound_analysis.items():
            activity = "shows" if props["is_chiral"] else "does not show"
            reason_str += f"- Compound {num} ({props['name']}) {activity} optical isomerism. Reason: {props['reason']}\n"
        
        reason_str += f"\nBased on this analysis, the compounds that show optical isomerism are {sorted(list(correct_compounds))}.\n"
        reason_str += f"This corresponds to option {correct_option}, which represents the set {sorted(list(options[correct_option]))}.\n"
        reason_str += f"The provided answer chose option {llm_option}, which represents the set {sorted(list(options[llm_option]))}, and is therefore incorrect."
        
        return reason_str

# The final answer provided by the user to be checked
llm_answer_to_check = """
To determine which of the given compounds show optical isomerism, we must analyze the structure of each molecule for chirality. A molecule is chiral, and therefore optically active, if it is non-superimposable on its mirror image. This typically arises from the presence of a chiral center or other elements of chirality, such as a chiral axis, and the absence of symmetry elements like a plane of symmetry or a center of inversion.

Here is a step-by-step analysis of each compound:

1.  **dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate**
    *   **Structure:** This is a biphenyl system where two phenyl rings are joined by a single C-C bond. All four positions adjacent to this bond (the *ortho* positions 2, 2', 6, and 6') are substituted with bulky groups (two nitro groups and two dimethyl carboxylate groups).
    *   **Careful Point (Atropisomerism):** In such substituted biphenyls, the large size of the ortho-substituents creates significant steric hindrance, which prevents the two rings from rotating freely around the central single bond. This restricted rotation locks the molecule into a twisted, non-planar conformation.
    *   **Analysis:** This locked conformation lacks both a plane of symmetry and a center of inversion. As a result, the molecule is non-superimposable on its mirror image. This phenomenon is known as atropisomerism, a form of axial chirality.
    *   **Conclusion:** The compound is chiral and exhibits optical isomerism.

2.  **methyl 2-hydroxypropanoate**
    *   **Structure:** Also known as methyl lactate, its structure is CH₃-CH(OH)-COOCH₃.
    *   **Careful Point (Chiral Center):** We must check for carbon atoms bonded to four different groups.
    *   **Analysis:** The carbon atom at position 2 is bonded to four distinct groups:
        1.  a hydrogen atom (-H)
        2.  a hydroxyl group (-OH)
        3.  a methyl group (-CH₃)
        4.  a methoxycarbonyl group (-COOCH₃)
    *   Since this carbon has four different substituents, it is a chiral center. A molecule containing a single chiral center is inherently chiral.
    *   **Conclusion:** The compound is chiral and exhibits optical isomerism.

3.  **benzophenone**
    *   **Structure:** Its chemical formula is (C₆H₅)₂C=O. It consists of a central carbonyl group attached to two phenyl rings.
    *   **Careful Point (Symmetry):** We need to check for elements of symmetry.
    *   **Analysis:** The molecule does not have a chiral center. It possesses a plane of symmetry that passes through the carbon-oxygen double bond and bisects the angle between the two phenyl rings. Any molecule with a plane of symmetry is achiral.
    *   **Conclusion:** The compound is achiral and does not exhibit optical isomerism.

4.  **dimethyl fumarate**
    *   **Structure:** This is the dimethyl ester of fumaric acid, the *trans*-isomer of but-2-enedioic acid.
    *   **Careful Point (Symmetry in Geometric Isomers):** We must check for symmetry elements in this *trans*-alkene.
    *   **Analysis:** The molecule is planar and possesses a center of inversion (a point of symmetry) at the midpoint of the carbon-carbon double bond. The presence of a center of inversion makes a molecule achiral.
    *   **Conclusion:** The compound is achiral and does not exhibit optical isomerism.

**Final Decision:**
Based on the analysis, only compounds **1** (dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate) and **2** (methyl 2-hydroxypropanoate) are chiral and therefore show optical isomerism. This corresponds to option A.

<<<A>>>
"""

# Run the check
result = check_chemistry_answer(llm_answer_to_check)
print(result)