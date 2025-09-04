import re

def check_optical_activity():
    """
    This function checks the correctness of the provided LLM answer regarding the number of optically active compounds.
    It does this by defining the ground truth for each compound based on chemical principles and comparing it to the LLM's conclusion.
    """

    # Ground truth analysis of each compound based on stereochemistry principles.
    # A compound is optically active if it is chiral and not a meso compound.
    compounds_analysis = [
        {
            "name": "(Z)-1-chloro-2-methylbut-1-ene",
            "is_active": False,
            "reason": "The molecule is planar around the C=C double bond, creating a plane of symmetry. It is achiral."
        },
        {
            "name": "(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione",
            "is_active": True,
            "reason": "The name specifies a single enantiomer (using R/S descriptors) of a complex, asymmetric molecule. It is chiral."
        },
        {
            "name": "(2R,3S)-2,3-dimethylsuccinic acid",
            "is_active": False,
            "reason": "This is a meso compound. Despite having two chiral centers, it has an internal element of symmetry (plane or center) and is achiral."
        },
        {
            "name": "(2R,3R)-2,3-dimethylsuccinic acid",
            "is_active": True,
            "reason": "This is a chiral diastereomer. It lacks internal symmetry and is non-superimposable on its mirror image (the S,S form)."
        },
        {
            "name": "(R)-cyclohex-3-en-1-ol",
            "is_active": True,
            "reason": "The name specifies a single enantiomer (R-form) of a chiral molecule. The carbon at position 1 is a chiral center."
        },
        {
            "name": "(1s,3s,5s)-cyclohexane-1,3,5-triol",
            "is_active": False,
            "reason": "This is the all-cis isomer, which is highly symmetric and possesses multiple planes of symmetry. It is achiral."
        },
        {
            "name": "1-cyclopentyl-3-methylbutan-1-one",
            "is_active": False,
            "reason": "The molecule has no chiral centers. The carbon at position 3 of the butane chain is bonded to two identical methyl groups. It is achiral."
        }
    ]

    # The correct count of optically active compounds based on our analysis.
    correct_count = sum(1 for c in compounds_analysis if c["is_active"])
    correct_active_compounds = {c["name"] for c in compounds_analysis if c["is_active"]}

    # The answer provided by the LLM.
    llm_answer_text = """
    To determine the number of optically active compounds, we must analyze each molecule for chirality. A compound is optically active if it is chiral and not present as a racemic mixture. A molecule is chiral if it is non-superimposable on its mirror image. The key is to look for chiral centers and the absence of symmetry elements like a plane of symmetry or a center of inversion.

    **Careful Points to Consider:**
    *   **Chirality:** A molecule is chiral if it lacks any improper axis of rotation (Sₙ), which includes planes of symmetry (σ = S₁) and centers of inversion (i = S₂).
    *   **Chiral Center:** A carbon atom bonded to four different groups is a common, but not exclusive, source of chirality.
    *   **Meso Compounds:** A molecule with chiral centers can still be achiral if it possesses an internal element of symmetry. Such a compound is called meso and is optically inactive.
    *   **Stereochemical Descriptors:** IUPAC names that include descriptors like (R), (S), or (E) specify a particular stereoisomer. If the name describes a single enantiomer of a chiral compound, it is optically active.

    ---

    ### **Step-by-Step Analysis of Each Compound**

    1.  **(Z)-1-chloro-2-methylbut-1-ene**
        *   **Analysis:** This molecule has a C=C double bond. The atoms involved in and directly attached to the double bond lie in a single plane. This molecular plane acts as a plane of symmetry. The molecule has no chiral centers.
        *   **Conclusion:** Achiral and **optically inactive**.

    2.  **(3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione**
        *   **Analysis:** The name includes the stereochemical descriptors `(3aR,7aS)`, which specify the absolute configuration at two chiral centers. This indicates that a single, specific stereoisomer is being described. The complex, bridged structure is asymmetric and lacks internal symmetry. A single enantiomer of a chiral compound is by definition optically active.
        *   **Conclusion:** Chiral and **optically active**.

    3.  **(2R,3S)-2,3-dimethylsuccinic acid**
        *   **Analysis:** This molecule has two chiral centers (C2 and C3) and is symmetric, with identical `HOOC-CH(CH3)-` groups. The `(2R,3S)` configuration results in a molecule with an internal element of symmetry (a center of inversion in the staggered conformation and a plane of symmetry in the eclipsed conformation). This is a classic example of a **meso compound**.
        *   **Conclusion:** Achiral (meso) and **optically inactive**.

    4.  **(2R,3R)-2,3-dimethylsuccinic acid**
        *   **Analysis:** This is a diastereomer of the meso compound. The `(2R,3R)` configuration lacks any internal symmetry elements. It is a chiral molecule, and its non-superimposable mirror image is the `(2S,3S)` isomer. The name specifies a single enantiomer.
        *   **Conclusion:** Chiral and **optically active**.

    5.  **(R)-cyclohex-3-en-1-ol**
        *   **Analysis:** The name specifies the `(R)` configuration at C1. The carbon at C1 is a chiral center because it is bonded to four different groups: -H, -OH, and two different paths around the ring (`-CH=CH-CH₂-` and `-CH₂-CH₂-`). The name refers to a single enantiomer of a chiral molecule.
        *   **Conclusion:** Chiral and **optically active**.

    6.  **(1s,3s,5s)-cyclohexane-1,3,5-triol**
        *   **Analysis:** This name describes the all-cis isomer, where all three hydroxyl groups are on the same side of the ring. This molecule is highly symmetric. It possesses three planes of symmetry, each passing through one C-OH group and the opposite CH₂ group. Due to these symmetry planes, the molecule is achiral.
        *   **Conclusion:** Achiral and **optically inactive**.

    7.  **1-cyclopentyl-3-methylbutan-1-one**
        *   **Analysis:** The structure is `(CH₃)₂CH-CH₂-C(=O)-cyclopentyl`. A search for chiral centers reveals none. The carbon at position 3 of the butane chain is bonded to two identical methyl groups, so it is not a chiral center. The molecule is achiral.
        *   **Conclusion:** Achiral and **optically inactive**.

    ### **Final Count**

    The compounds that will exhibit optical activity are:
    1.  (3aR,7aS,E)-8-(chloromethylene)hexahydro-4,7-methanoisobenzofuran-1,3-dione
    2.  (2R,3R)-2,3-dimethylsuccinic acid
    3.  (R)-cyclohex-3-en-1-ol

    There are a total of **3** optically active compounds. This corresponds to option B.

    <<<B>>>
    """

    # Extract the final letter answer and the count from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find the final answer in the format <<<X>>>."

    llm_letter_answer = match.group(1)
    answer_map = {'A': 5, 'B': 3, 'C': 4, 'D': 2}
    llm_count = answer_map.get(llm_letter_answer)

    # Check if the count is correct
    if llm_count != correct_count:
        return f"Incorrect. The final answer is {llm_count}, but the correct number of optically active compounds is {correct_count}."

    # Check if the specific compounds identified as active are correct
    llm_active_compounds_text = re.search(r'The compounds that will exhibit optical activity are:(.*?)There are a total of', llm_answer_text, re.DOTALL)
    if not llm_active_compounds_text:
        return "Could not parse the list of active compounds from the answer text for verification."
    
    llm_active_list = llm_active_compounds_text.group(1).strip().split('\n')
    # Clean up the list extracted from the text
    llm_identified_active = set()
    for item in llm_active_list:
        # Find the first compound name in the line
        for compound in compounds_analysis:
            if compound["name"] in item:
                llm_identified_active.add(compound["name"])
                break

    if llm_identified_active != correct_active_compounds:
        missing = correct_active_compounds - llm_identified_active
        extra = llm_identified_active - correct_active_compounds
        error_msg = "Incorrect. The count is right, but the specific compounds identified are wrong.\n"
        if missing:
            error_msg += f"Missing from LLM's list: {', '.join(missing)}\n"
        if extra:
            error_msg += f"Incorrectly included in LLM's list: {', '.join(extra)}\n"
        return error_msg

    # If both the count and the specific compounds are correct
    return "Correct"

# Run the check and print the result
result = check_optical_activity()
print(result)