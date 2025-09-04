import re

def check_chirality(compound_name):
    """
    Determines if a given compound is chiral based on its name and known chemical rules.
    Returns a tuple: (is_chiral, reason).
    """
    if compound_name == "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate":
        # This is a classic example of atropisomerism. Bulky groups at all four ortho positions
        # (2, 2', 6, 6') restrict rotation around the central C-C bond, creating a stable,
        # non-planar, chiral molecule that lacks a plane of symmetry.
        return (True, "Exhibits optical isomerism due to atropisomerism (axial chirality).")
    
    elif compound_name == "methyl 2-hydroxypropanoate":
        # The structure is CH3-CH(OH)-COOCH3. The carbon at position 2 is bonded to four
        # different groups: -H, -OH, -CH3, and -COOCH3. This makes it a chiral center.
        return (True, "Exhibits optical isomerism due to a chiral center.")
        
    elif compound_name == "benzophenone":
        # The structure is (C6H5)2C=O. The molecule has a plane of symmetry that passes
        # through the C=O bond and bisects the angle between the two phenyl rings.
        # Molecules with a plane of symmetry are achiral.
        return (False, "Does not exhibit optical isomerism because it is achiral (has a plane of symmetry).")
        
    elif compound_name == "dimethyl fumarate":
        # This is the trans-isomer of an alkene. The molecule is planar and has a center of
        # inversion at the midpoint of the C=C double bond. Both features make it achiral.
        return (False, "Does not exhibit optical isomerism because it is achiral (planar and has a center of inversion).")
        
    else:
        return (None, "Unknown compound.")

def check_correctness_of_answer(llm_answer):
    """
    Checks if the final answer from an LLM is correct for the given chemistry question.

    Args:
        llm_answer: The string output from the LLM, which should contain the final answer
                    in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a reason for the incorrectness.
    """
    # Extract the final answer choice from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Failure: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."
    
    answer_choice = match.group(1)

    # Define the compounds from the question
    compounds = [
        "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
        "methyl 2-hydroxypropanoate",
        "benzophenone",
        "dimethyl fumarate"
    ]

    # Determine the correct set of optically active compounds based on chemical principles
    correct_optically_active_set = set()
    analysis_log = []
    for i, compound in enumerate(compounds, 1):
        is_chiral, reason = check_chirality(compound)
        analysis_log.append(f"Compound {i} ({compound}): {reason}")
        if is_chiral:
            correct_optically_active_set.add(i)

    # Define the sets corresponding to each multiple-choice option
    options = {
        'A': {2, 3},
        'B': {1, 2},
        'C': {1, 2, 4},
        'D': {3, 4}
    }

    # Get the set of compounds corresponding to the LLM's answer
    llm_answer_set = options.get(answer_choice)

    # Compare the LLM's answer set with the correct set
    if llm_answer_set == correct_optically_active_set:
        return "Correct"
    else:
        # Construct a detailed reason for the incorrectness
        reason = (f"Incorrect. The provided answer is '{answer_choice}', which corresponds to compounds {sorted(list(llm_answer_set))}.\n"
                  f"The correct answer is 'B', which corresponds to compounds {sorted(list(correct_optically_active_set))}.\n\n"
                  "Reasoning:\n")
        
        # Check for compounds incorrectly included
        incorrectly_included = llm_answer_set - correct_optically_active_set
        if incorrectly_included:
            for i in sorted(list(incorrectly_included)):
                _, analysis = check_chirality(compounds[i-1])
                reason += f"- The answer incorrectly includes compound {i}, which is not optically active. {analysis}\n"

        # Check for compounds incorrectly excluded
        incorrectly_excluded = correct_optically_active_set - llm_answer_set
        if incorrectly_excluded:
            for i in sorted(list(incorrectly_excluded)):
                _, analysis = check_chirality(compounds[i-1])
                reason += f"- The answer incorrectly excludes compound {i}, which is optically active. {analysis}\n"
                
        return reason.strip()

# The user has provided the final answer from the LLM in the prompt.
# We will use this as the input for our checking function.
final_llm_answer = """
Based on an analysis of the provided answers and a step-by-step evaluation of the chemical principles, here is the breakdown to determine the correct answer.

**Objective**: Identify which of the four given compounds exhibit optical isomerism. A compound shows optical isomerism if its molecular structure is chiral, meaning it is non-superimposable on its mirror image.

**Step-by-step Analysis of Each Compound:**

1.  **dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate**
    *   **Careful Point**: This is a substituted biphenyl. The key feature to check for is atropisomerism, a type of axial chirality that arises from restricted rotation around a single bond.
    *   **Analysis**: The molecule has bulky substituents (nitro -NO₂ and dimethyl carboxylate -COOCH₃) at all four *ortho* positions (2, 2', 6, and 6'). This significant steric hindrance prevents the two phenyl rings from rotating freely and becoming coplanar. The molecule is locked into a stable, twisted conformation. This conformation lacks a plane of symmetry or a center of inversion, making the molecule chiral.
    *   **Conclusion**: This compound exhibits optical isomerism.

2.  **methyl 2-hydroxypropanoate**
    *   **Careful Point**: The most common source of chirality is a chiral center, which is a carbon atom bonded to four different groups.
    *   **Analysis**: The structure is CH₃-CH(OH)-COOCH₃. The carbon atom at position 2 is bonded to four distinct groups:
        1.  -H (hydrogen)
        2.  -OH (hydroxyl)
        3.  -CH₃ (methyl)
        4.  -COOCH₃ (methoxycarbonyl)
    *   **Conclusion**: Since it contains a chiral center, the molecule is chiral and exhibits optical isomerism.

3.  **benzophenone**
    *   **Careful Point**: Molecules with elements of symmetry, such as a plane of symmetry or a center of inversion, are achiral.
    *   **Analysis**: The structure is (C₆H₅)₂C=O. The molecule possesses a plane of symmetry that passes through the central C=O bond and bisects the angle between the two identical phenyl rings.
    *   **Conclusion**: Due to the plane of symmetry, the molecule is achiral and does not show optical isomerism.

4.  **dimethyl fumarate**
    *   **Careful Point**: This is a geometric isomer (*trans*). We must check for chirality, not just stereoisomerism in general.
    *   **Analysis**: The structure is the *trans* isomer of an alkene. The molecule is planar and has a center of inversion at the midpoint of the carbon-carbon double bond. The presence of either a plane of symmetry (the molecular plane itself) or a center of inversion makes a molecule achiral.
    *   **Conclusion**: The molecule is achiral and does not show optical isomerism.

**Final Decision:**

*   Compound 1 is optically active.
*   Compound 2 is optically active.
*   Compound 3 is not optically active.
*   Compound 4 is not optically active.

Therefore, the correct answer must include compounds 1 and 2. Reviewing the options:
A) 2 and 3
B) 1 and 2
C) 1, 2 and 4
D) 3 and 4

The correct option is B. This aligns with the chemical reasoning presented in the majority of the candidate answers, even though many of them had typographical errors in their final letter choice.

<<<B>>>
"""

# Run the check
result = check_correctness_of_answer(final_llm_answer)
print(result)