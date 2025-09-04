import re

def check_chemistry_answer(answer_text):
    """
    Checks the correctness of the LLM's answer for the optical activity question.
    """
    # --- Ground Truth Definition ---
    # Based on chemical principles of chirality and symmetry.
    # True = Optically Active, False = Optically Inactive
    ground_truth_activity = {
        "1": False, # (Z)-1-chloro-2-methylbut-1-ene (achiral, plane of symmetry)
        "2": True,  # (3aR,7aS,E)-... (chiral, specified enantiomer)
        "3": False, # (2R,3S)-2,3-dimethylsuccinic acid (achiral, meso compound)
        "4": True,  # (2R,3R)-2,3-dimethylsuccinic acid (chiral, specified enantiomer)
        "5": True,  # (R)-cyclohex-3-en-1-ol (chiral, specified enantiomer)
        "6": False, # (1s,3s,5s)-cyclohexane-1,3,5-triol (achiral, planes of symmetry)
        "7": False  # 1-cyclopentyl-3-methylbutan-1-one (achiral, no chiral center)
    }
    correct_active_count = sum(1 for is_active in ground_truth_activity.values() if is_active)

    # Question options
    options_map = {
        "A": 4,
        "B": 2,
        "C": 5,
        "D": 3
    }

    # --- Parsing the LLM's Answer ---
    # 1. Extract the final selected option (e.g., <<<D>>>)
    final_answer_match = re.search(r'<<<([A-D])>>>', answer_text)
    if not final_answer_match:
        return "Constraint not satisfied: The answer does not contain a final answer in the format <<<X>>>."
    selected_option = final_answer_match.group(1)

    # 2. Extract the number of active compounds claimed in the text
    # Use a flexible regex to find phrases like "total of 3", "3 optically active", etc.
    count_match = re.search(r'total of \D*(\d+)', answer_text, re.IGNORECASE)
    if not count_match:
        return "Constraint not satisfied: Could not find the total count of active compounds stated in the answer's text."
    claimed_count = int(count_match.group(1))

    # --- Verification ---
    # Check 1: Does the claimed count in the text match the ground truth?
    if claimed_count != correct_active_count:
        return f"Incorrect: The answer claims there are {claimed_count} optically active compounds, but the correct number is {correct_active_count}."

    # Check 2: Does the selected option correspond to the correct count?
    value_of_selected_option = options_map.get(selected_option)
    if value_of_selected_option != correct_active_count:
        return f"Incorrect: The answer selected option {selected_option} which corresponds to {value_of_selected_option} compounds, but the correct number of active compounds is {correct_active_count}."

    # Check 3: Internal consistency check (does the claimed count match the selected option?)
    if claimed_count != value_of_selected_option:
        return f"Incorrect: The answer's text claims there are {claimed_count} active compounds, but it selects option {selected_option} which corresponds to {value_of_selected_option}."

    # If all checks pass, the answer is correct.
    return "Correct"

# The user-provided text containing the question and the final answer to be checked.
llm_answer_to_check = """
To determine the number of optically active compounds, we must analyze each molecule for chirality. A compound is optically active if it is chiral (non-superimposable on its mirror image) and is not present as a 50:50 racemic mixture. The presence of stereochemical descriptors like (R) or (S) in a name indicates that a single, specific stereoisomer is being considered.

**Careful Points for Analysis:**
*   **Chirality:** The primary requirement for optical activity. A molecule is chiral if it lacks any improper axis of rotation (Sₙ), which includes planes of symmetry (σ = S₁) and centers of inversion (i = S₂).
*   **Chiral Centers:** A carbon atom bonded to four different groups is a common source of chirality, but its presence does not guarantee the entire molecule is chiral.
*   **Meso Compounds:** These are molecules that contain chiral centers but are achiral overall due to an internal element of symmetry. They are optically inactive.
*   **Symmetry:** A thorough search for planes of symmetry and centers of inversion is crucial, especially in cyclic and seemingly complex molecules.

---

### **Step-by-Step Analysis of Each Compound**

1.  **(Z)-1-chloro-2-methylbut-1-ene**:
    *   **Analysis**: This is an alkene. The atoms involved in the double bond and the atoms directly attached to them lie in a single plane. This molecular plane acts as a plane of symmetry.
    *   **Conclusion**: Achiral and **optically inactive**.

2.  **(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione**:
    *   **Analysis**: The name includes the stereochemical descriptors `(3aR,7aS)`, which specify the absolute configuration at two chiral centers. This indicates a single enantiomer of a complex, chiral molecule is being described. The bridged structure is asymmetric.
    *   **Conclusion**: Chiral and **optically active**.

3.  **(2R,3S)-2,3-dimethylsuccinic acid**:
    *   **Analysis**: This molecule has two chiral centers (C2 and C3) but is symmetric, with identical groups at both ends (`HOOC-CH(CH3)-`). The `(2R,3S)` configuration results in an internal element of symmetry (a plane or center of inversion), making it a meso compound.
    *   **Conclusion**: Achiral (meso) and **optically inactive**.

4.  **(2R,3R)-2,3-dimethylsuccinic acid**:
    *   **Analysis**: This is a diastereomer of the meso compound. The `(2R,3R)` configuration lacks any internal symmetry elements. Its mirror image is the non-superimposable `(2S,3S)` isomer. The name specifies a single enantiomer.
    *   **Conclusion**: Chiral and **optically active**.

5.  **(R)-cyclohex-3-en-1-ol**:
    *   **Analysis**: The name specifies the `(R)` configuration at C1. This carbon is a chiral center because it is bonded to four different groups: -H, -OH, and two different paths around the ring (`-CH=CH-CH₂-` and `-CH₂-CH₂-CH=`). The name refers to a single enantiomer.
    *   **Conclusion**: Chiral and **optically active**.

6.  **(1s,3s,5s)-cyclohexane-1,3,5-triol**:
    *   **Analysis**: This notation describes the all-cis isomer, where all three hydroxyl groups are on the same side of the ring. This molecule is highly symmetric and possesses three planes of symmetry, each passing through one C-OH group and the opposite CH₂ group.
    *   **Conclusion**: Achiral and **optically inactive**.

7.  **1-cyclopentyl-3-methylbutan-1-one**:
    *   **Analysis**: A search for chiral centers reveals none. The carbon at position 3 of the butane chain is bonded to two identical methyl groups, so it is not a chiral center. The molecule as a whole is achiral.
    *   **Conclusion**: Achiral and **optically inactive**.

### **Final Count**

The compounds that will exhibit optical activity are:
*   (3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione
*   (2R,3R)-2,3-dimethylsuccinic acid
*   (R)-cyclohex-3-en-1-ol

There are a total of **3** optically active compounds. This corresponds to option D.

<<<D>>>
"""

# Run the check
result = check_chemistry_answer(llm_answer_to_check)
print(result)