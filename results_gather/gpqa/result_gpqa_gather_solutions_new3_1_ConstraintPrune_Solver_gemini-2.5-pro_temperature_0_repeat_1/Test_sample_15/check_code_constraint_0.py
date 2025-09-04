import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by comparing it to a ground truth
    analysis of each compound's optical activity.
    """
    # Ground truth based on established chemical principles.
    # A compound is optically active if it is chiral and not a meso compound.
    ground_truth = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "is_active": False,
            "reason": "Achiral. Has a plane of symmetry (the plane of the double bond). No chiral centers."
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "is_active": True,
            "reason": "Chiral. The name specifies a single enantiomer with defined stereocenters (3aR, 7aS)."
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "is_active": False,
            "reason": "Achiral (meso compound). Has two chiral centers but also an internal plane of symmetry."
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "is_active": True,
            "reason": "Chiral. This is a single enantiomer of a chiral pair ((2R,3R) and (2S,3S)). It lacks internal symmetry."
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "is_active": True,
            "reason": "Chiral. The name specifies a single enantiomer (R) at a chiral center (C1)."
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "is_active": False,
            "reason": "Achiral. The all-cis configuration results in a highly symmetric molecule with multiple planes of symmetry."
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "is_active": False,
            "reason": "Achiral. No chiral centers in the molecule."
        }
    ]

    # The final answer provided by the LLM.
    llm_answer_text = """
    To decide the final answer, I will analyze each compound for optical activity. A compound is optically active if it is chiral and not present as a 50:50 mixture of enantiomers (a racemic mixture). A molecule is chiral if it is non-superimposable on its mirror image. Key features to look for are chiral centers and elements of symmetry (like a plane of symmetry or a center of inversion).

    **Step-by-step analysis of each compound:**

    1.  **(Z)-1-chloro-2-methylbut-1-ene**:
        *   **Structure**: `Cl-CH=C(CH3)(CH2CH3)`
        *   **Analysis**: This molecule is an alkene. The atoms involved in the double bond and the atoms directly attached to them lie in a single plane. This plane acts as a plane of symmetry for the entire molecule. The molecule has no chiral centers.
        *   **Conclusion**: Achiral and **optically inactive**.

    2.  **(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione**:
        *   **Analysis**: The name includes specific stereochemical descriptors: `(3aR,7aS,E)`. This indicates that a single, specific stereoisomer is being described. The presence of `R` and `S` confirms the existence of chiral centers. A complex molecule like this, specified as a single stereoisomer, is chiral.
        *   **Conclusion**: Chiral and **optically active**.

    3.  **(2R,3S)-2,3-dimethylsuccinic acid**:
        *   **Structure**: `HOOC-CH(CH3)-CH(CH3)-COOH`
        *   **Analysis**: This molecule has two chiral centers (C2 and C3). However, the molecule is symmetric about its midpoint. The `(2R,3S)` configuration results in an internal plane of symmetry. A molecule with chiral centers but also a plane of symmetry is called a meso compound.
        *   **Conclusion**: Achiral (meso) and **optically inactive**.

    4.  **(2R,3R)-2,3-dimethylsuccinic acid**:
        *   **Analysis**: This is a diastereomer of the meso compound above. The `(2R,3R)` configuration does not have an internal plane of symmetry. Its mirror image is the `(2S,3S)` isomer, and they are non-superimposable. It is a chiral molecule. Since a single enantiomer is specified, it is optically active.
        *   **Conclusion**: Chiral and **optically active**.

    5.  **(R)-cyclohex-3-en-1-ol**:
        *   **Analysis**: The name specifies the `(R)` configuration at a chiral center (C1, the carbon with the -OH group). This carbon is bonded to four different groups: -H, -OH, and two different pathways around the ring. The molecule is chiral, and a single enantiomer is specified.
        *   **Conclusion**: Chiral and **optically active**.

    6.  **(1s,3s,5s)-cyclohexane-1,3,5-triol**:
        *   **Analysis**: This name describes the all-cis isomer, where all three hydroxyl groups are on the same side of the ring. This molecule is highly symmetric. It possesses three planes of symmetry, each passing through one C-OH group and the opposite CH2 group.
        *   **Conclusion**: Achiral and **optically inactive**.

    7.  **1-cyclopentyl-3-methylbutan-1-one**:
        *   **Structure**: `(CH3)2CH-CH2-C(=O)-Cyclopentyl`
        *   **Analysis**: We check for chiral centers. The carbon at position 3 is bonded to a hydrogen, two identical methyl groups, and a `-CH2-C(=O)-cyclopentyl` group. Since it is bonded to two identical groups, it is not a chiral center. The molecule lacks any chiral centers and is achiral.
        *   **Conclusion**: Achiral and **optically inactive**.

    **Final Count:**
    The compounds that will exhibit optical activity are:
    1.  (3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione
    2.  (2R,3R)-2,3-dimethylsuccinic acid
    3.  (R)-cyclohex-3-en-1-ol

    There are a total of **3** optically active compounds.

    Matching this count to the options:
    A) 2
    B) 4
    C) 3
    D) 5

    The correct option is C.

    <<<C>>>
    """

    # 1. Calculate the correct count from ground truth
    correct_count = sum(1 for c in ground_truth if c["is_active"])
    correct_active_compounds = {c["name"] for c in ground_truth if c["is_active"]}

    # 2. Parse the LLM's final answer
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not find the final answer in the format <<<...>>> in the provided text."
    
    llm_option = match.group(1)
    option_map = {'A': 2, 'B': 4, 'C': 3, 'D': 5}
    llm_count = option_map.get(llm_option)

    # 3. Compare the counts
    if llm_count != correct_count:
        return (f"Incorrect. The answer states there are {llm_count} optically active compounds, "
                f"but the correct number is {correct_count}.")

    # 4. Check the reasoning by comparing the list of active compounds
    # Extract the list of active compounds from the LLM's reasoning
    reasoning_match = re.search(r'The compounds that will exhibit optical activity are:\s*\n(.*?)\n\n', llm_answer_text, re.DOTALL)
    if not reasoning_match:
        return "Failure: Could not parse the list of active compounds from the reasoning text."
    
    llm_active_list_text = reasoning_match.group(1).strip()
    # Split by newline and remove the leading number like "1. "
    llm_active_compounds = {re.sub(r'^\d+\.\s*', '', line.strip()) for line in llm_active_list_text.split('\n')}

    if llm_active_compounds != correct_active_compounds:
        missing = correct_active_compounds - llm_active_compounds
        extra = llm_active_compounds - correct_active_compounds
        error_msg = "Incorrect. The final count is correct, but the reasoning is flawed.\n"
        if missing:
            error_msg += f"The reasoning incorrectly omits the following active compounds: {list(missing)}\n"
        if extra:
            error_msg += f"The reasoning incorrectly includes the following inactive compounds as active: {list(extra)}\n"
        return error_msg.strip()

    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)