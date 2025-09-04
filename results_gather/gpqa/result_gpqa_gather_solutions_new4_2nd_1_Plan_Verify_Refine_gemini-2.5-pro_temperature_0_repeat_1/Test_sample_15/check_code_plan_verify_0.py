import re

def check_correctness(answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the optical activity question.

    The function verifies:
    1. The analysis of each individual compound against established chemical principles.
    2. The final count of optically active compounds.
    3. The consistency between the final count and the chosen option letter.
    """

    # --- Ground Truth ---
    # This dictionary represents the correct chemical analysis for each compound.
    # 'active': True if optically active, False otherwise.
    # 'reason': A brief explanation for the classification.
    ground_truth = {
        "(Z)-1-chloro-2-methylbut-1-ene": {
            "active": False,
            "reason": "Achiral due to a plane of symmetry (planar alkene)."
        },
        "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione": {
            "active": True,
            "reason": "Chiral; the name specifies a single enantiomer of a complex, asymmetric molecule."
        },
        "(2R,3S)-2,3-dimethylsuccinic acid": {
            "active": False,
            "reason": "Achiral meso compound with an internal plane of symmetry."
        },
        "(2R,3R)-2,3-dimethylsuccinic acid": {
            "active": True,
            "reason": "Chiral; the name specifies a single enantiomer (diastereomer of the meso form)."
        },
        "(R)-cyclohex-3-en-1-ol": {
            "active": True,
            "reason": "Chiral; the name specifies a single enantiomer with a chiral center at C1."
        },
        "(1s,3s,5s)-cyclohexane-1,3,5-triol": {
            "active": False,
            "reason": "Achiral due to multiple planes of symmetry (all-cis isomer)."
        },
        "1-cyclopentyl-3-methylbutan-1-one": {
            "active": False,
            "reason": "Achiral; lacks any chiral centers."
        }
    }

    # --- Parse the LLM's Answer ---
    # Extract the LLM's conclusion for each compound.
    try:
        # Split the answer into sections, one for each compound analysis
        sections = re.split(r'\n\s*\d\.\s+', answer_text)
        answer_analysis = {}
        
        for compound_name in ground_truth.keys():
            # Escape special characters in the compound name for safe regex searching
            safe_name = re.escape(compound_name)
            found_compound = False
            for section in sections:
                if re.search(safe_name, section):
                    found_compound = True
                    # Find the conclusion line and determine if it says "active" or "inactive"
                    conclusion_match = re.search(r'Conclusion\s*:\s*.*?optically\s+(active|inactive)', section, re.DOTALL | re.IGNORECASE)
                    if conclusion_match:
                        status = conclusion_match.group(1).lower()
                        answer_analysis[compound_name] = (status == "active")
                    else:
                        return f"Constraint not satisfied: Could not parse the conclusion for '{compound_name}'."
                    break
            if not found_compound:
                 return f"Constraint not satisfied: The answer did not provide an analysis for the compound '{compound_name}'."

    except Exception as e:
        return f"Failed to parse the answer's analysis for each compound. Error: {e}"

    # --- Verification Steps ---

    # 1. Verify the analysis of each compound
    for compound, data in ground_truth.items():
        if answer_analysis[compound] != data["active"]:
            is_or_is_not_answer = "is" if answer_analysis[compound] else "is not"
            is_or_is_not_truth = "is" if data["active"] else "is not"
            return f"Constraint not satisfied: Incorrect analysis for '{compound}'. The answer claims it {is_or_is_not_answer} optically active, but it {is_or_is_not_truth}. Reason: {data['reason']}"

    # 2. Verify the final count
    correct_count = sum(1 for data in ground_truth.values() if data["active"])
    
    count_match = re.search(r"total of \*\*(\d+)\*\* optically active compounds", answer_text)
    if not count_match:
        # Fallback for answers without markdown
        count_match = re.search(r"total of (\d+) optically active compounds", answer_text)
    
    if not count_match:
        return "Constraint not satisfied: Could not find the final count of active compounds in the answer text."
    
    answer_count = int(count_match.group(1))
    
    if answer_count != correct_count:
        return f"Constraint not satisfied: The final count is incorrect. The answer states there are {answer_count} optically active compounds, but the correct number is {correct_count}."

    # 3. Verify the final letter choice and its consistency with the options
    options = {'A': 3, 'B': 5, 'C': 4, 'D': 2}
    
    final_answer_match = re.search(r"<<<([A-D])>>>", answer_text)
    if not final_answer_match:
        return "Constraint not satisfied: Could not find the final answer choice in the format <<<X>>>."
        
    final_answer_letter = final_answer_match.group(1)
    
    if options.get(final_answer_letter) != answer_count:
        return f"Constraint not satisfied: The final answer choice '{final_answer_letter}' is inconsistent with the answer's own derived count of {answer_count}. Option {final_answer_letter} corresponds to {options.get(final_answer_letter)} compounds."

    # If all checks pass, the answer is correct.
    return "Correct"

# The final answer provided by the user to be checked
llm_answer_to_check = """
To determine the final answer, a systematic analysis of each compound for chirality is required. A compound exhibits optical activity if it is chiral (non-superimposable on its mirror image) and is present as a single enantiomer (not a racemic mixture).

### **Careful Points for Analysis**

*   **Chirality:** The fundamental requirement for optical activity. A molecule is chiral if it lacks any improper axis of rotation (Sₙ), which includes common symmetry elements like a plane of symmetry (σ = S₁) and a center of inversion (i = S₂).
*   **Chiral Center:** A carbon atom bonded to four different groups is a common source of chirality, but its presence does not guarantee the entire molecule is chiral.
*   **Meso Compounds:** These are molecules that contain chiral centers but are achiral overall due to an internal element of symmetry. They are optically inactive.
*   **Stereochemical Descriptors:** IUPAC names that include descriptors like (R), (S), or specific locants like (3aR, 7aS) refer to a single, specific stereoisomer. If that isomer is chiral, it is considered optically active for the purpose of this question.

### **Step-by-Step Analysis of Each Compound**

1.  **(Z)-1-chloro-2-methylbut-1-ene**: This is a planar alkene. The plane of the double bond acts as a plane of symmetry for the entire molecule.
    *   **Conclusion**: Achiral and **optically inactive**.

2.  **(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione**: The name explicitly defines a single stereoisomer with specific configurations `(3aR, 7aS)` at its chiral centers. This complex, bridged molecule is asymmetric. As a single enantiomer of a chiral compound, it is optically active.
    *   **Conclusion**: Chiral and **optically active**.

3.  **(2R,3S)-2,3-dimethylsuccinic acid**: This molecule has two chiral centers (C2 and C3) but is structurally symmetric. The `(2R,3S)` configuration results in an internal plane of symmetry, making it a **meso compound**.
    *   **Conclusion**: Achiral (meso) and **optically inactive**.

4.  **(2R,3R)-2,3-dimethylsuccinic acid**: This is a diastereomer of the meso compound. The `(2R,3R)` configuration lacks any internal symmetry. It is a chiral molecule, and the name specifies a single enantiomer.
    *   **Conclusion**: Chiral and **optically active**.

5.  **(R)-cyclohex-3-en-1-ol**: The name specifies the `(R)` configuration at C1. This carbon is a chiral center because it is bonded to four different groups: -H, -OH, and two different paths around the ring. The name refers to a single enantiomer of a chiral molecule.
    *   **Conclusion**: Chiral and **optically active**.

6.  **(1s,3s,5s)-cyclohexane-1,3,5-triol**: This describes the all-cis isomer, where all three hydroxyl groups are on the same side of the ring. This highly symmetric molecule possesses three planes of symmetry.
    *   **Conclusion**: Achiral and **optically inactive**.

7.  **1-cyclopentyl-3-methylbutan-1-one**: The structure is `(CH₃)₂CH-CH₂-C(=O)-cyclopentyl`. There are no chiral centers in this molecule. The carbon at position 3 of the butane chain is bonded to two identical methyl groups.
    *   **Conclusion**: Achiral and **optically inactive**.

### **Final Count**

The compounds that will exhibit optical activity are:
1.  (3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione
2.  (2R,3R)-2,3-dimethylsuccinic acid
3.  (R)-cyclohex-3-en-1-ol

There are a total of **3** optically active compounds. According to the options provided in the question:
A) 3
B) 5
C) 4
D) 2

The correct option is A.

<<<A>>>
"""

# Run the check and print the result
result = check_correctness(llm_answer_to_check)
print(result)