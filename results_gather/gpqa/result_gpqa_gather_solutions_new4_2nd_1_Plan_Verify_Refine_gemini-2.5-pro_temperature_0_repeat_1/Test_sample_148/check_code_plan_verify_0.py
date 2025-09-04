import re

def check_chemistry_answer(llm_response_text: str) -> str:
    """
    Checks the correctness of the answer to the peptide analysis question.

    The function codifies the experimental observations and the properties of each
    possible explanation (precursor, enantiomers, diastereoisomers, etc.) to
    logically deduce the correct answer. It then compares this correct answer
    to the one provided by the LLM.

    Args:
        llm_response_text: The full text of the LLM's response, including the
                           final answer in the format <<<X>>>.

    Returns:
        "Correct" if the LLM's answer is correct.
        A string explaining the error if the LLM's answer is incorrect.
    """

    # 1. Define the key experimental observations from the question text.
    # These are the constraints that the correct answer must satisfy.
    observations = {
        "has_same_mass": True,  # From MS: "Both peaks have the same mass spectrum"
        "is_separable_by_lc": True, # From LC: "two clearly defined peaks"
        "is_distinguishable_by_nmr": True # From NMR: "two peaks that both correspond to the same alpha-proton"
    }

    # 2. Define the properties of each possible explanation based on chemical principles.
    # The analysis is assumed to be done in a standard (achiral) environment.
    explanations = {
        "A": {
            "description": "'Double coupling' product",
            "has_same_mass": False,  # A double-coupled product would have a higher mass.
            "is_separable_by_lc": True,
            "is_distinguishable_by_nmr": True
        },
        "B": {
            "description": "Contamination with a precursor",
            "has_same_mass": False,  # A precursor would have a different (likely lower) mass.
            "is_separable_by_lc": True,
            "is_distinguishable_by_nmr": True
        },
        "C": {
            "description": "Mixture of enantiomers",
            "has_same_mass": True,
            "is_separable_by_lc": False,  # Enantiomers are not separable by standard (achiral) LC.
            "is_distinguishable_by_nmr": False # Enantiomers are not distinguishable by standard (achiral) NMR.
        },
        "D": {
            "description": "Mixture of diastereoisomers",
            "has_same_mass": True,
            "is_separable_by_lc": True,
            "is_distinguishable_by_nmr": True
        }
    }

    # 3. Determine the logically correct option by finding which explanation matches all observations.
    correct_option = None
    for option_key, properties in explanations.items():
        if all(properties[obs_key] == obs_value for obs_key, obs_value in observations.items()):
            correct_option = option_key
            break

    # 4. Extract the final answer provided by the LLM.
    match = re.search(r'<<<([A-D])>>>', llm_response_text)
    if not match:
        return "Error: Could not find the final answer in the required format '<<<X>>>' in the provided text."
    
    llm_choice = match.group(1)

    # 5. Compare the LLM's answer with the logically derived correct answer.
    if llm_choice == correct_option:
        return "Correct"
    else:
        # If incorrect, provide a detailed explanation of the error.
        llm_props = explanations[llm_choice]
        correct_props = explanations[correct_option]
        
        reason = (f"Incorrect. The provided answer was '{llm_choice}', but the correct answer is '{correct_option}'.\n\n"
                  f"Here is the logical breakdown:\n"
                  f"1. **MS Data (Same Mass):** The two species have the same mass, meaning they are isomers. This rules out options A and B, which would have different masses.\n"
                  f"2. **LC & NMR Data (Separable/Distinguishable):** The two species are separated by LC and distinguished by NMR. This means they have different physical and chemical properties.\n"
                  f"3. **Conclusion:** Enantiomers (Option C) have identical properties in an achiral environment and would not be separated or distinguished. Diastereoisomers (Option D) are isomers with different properties, fitting all observations.\n\n"
                  f"The chosen answer '{llm_choice}' ({llm_props['description']}) is wrong because it fails the following checks:\n")
        
        if llm_props['has_same_mass'] != observations['has_same_mass']:
            reason += "- It is inconsistent with the MS observation that both species have the same mass.\n"
        if llm_props['is_separable_by_lc'] != observations['is_separable_by_lc']:
            reason += "- It is inconsistent with the LC observation that the species are separable.\n"
        if llm_props['is_distinguishable_by_nmr'] != observations['is_distinguishable_by_nmr']:
            reason += "- It is inconsistent with the NMR observation that the species are distinguishable.\n"
            
        return reason

# The user-provided response to be checked
llm_answer_text = """
My plan is to first deconstruct the analytical data provided in the question to establish a set of core facts. Then, I will systematically evaluate each of the four options (A, B, C, D) against these facts. I will also analyze the reasoning provided by the 17 candidate answers to find the logical consensus, noting why their final letter choices might differ. Finally, I will synthesize this information to arrive at the definitive correct answer.

**1. Deconstruction of the Experimental Data**

*   **¹H NMR Data**: The observation of two distinct peaks for a single alpha-proton, with equal integrals summing to 1H, is a key finding. Since spin-spin coupling is ruled out, this means the proton exists in two different chemical environments. The equal integrals indicate that the two molecular species responsible are present in a roughly 1:1 ratio.
*   **LC-MS Data**:
    *   **LC (Liquid Chromatography)**: The separation of the sample into two distinct peaks confirms the presence of two different chemical species that have different physical properties (e.g., polarity), allowing them to be separated.
    *   **MS (Mass Spectrometry)**: The fact that both peaks have the *same mass spectrum*, which is consistent with the expected product, is the most critical piece of evidence. It proves that the two species are **isomers**—they have the same chemical formula and mass but a different arrangement of atoms.

**2. Analysis of Candidate Answers**

A review of the 17 candidate answers reveals a powerful consensus in their *reasoning*, even though their final letter choices are scattered (A, B, C, and D are all selected by various agents). This discrepancy in the final letter strongly suggests that the options were presented in a different, randomized order to each agent. However, the logical pathway to the conclusion is identical in every detailed response:
*   All agents correctly use the MS data to conclude the species are isomers, thus ruling out precursors or side-products with different masses.
*   All agents correctly use the LC and NMR data (two separable peaks) to rule out enantiomers, which are indistinguishable in standard achiral environments.
*   All agents correctly conclude that the only remaining possibility consistent with all data is a mixture of diastereoisomers.

**3. Step-by-Step Evaluation of the Original Options**

Based on the established facts and the consensus reasoning, let's evaluate the original options:

*   **A) 'Double coupling' has occurred during an amide-bond forming reaction:** This would result in a product with a higher molecular weight. This is contradicted by the MS data showing both species have the expected mass. This is incorrect.
*   **B) The compound is contaminated with a precursor:** A precursor would have a different, likely lower, molecular weight. This is also contradicted by the MS data. This is incorrect.
*   **C) The crude compound exists as a mixture of enantiomers:** Enantiomers are non-superimposable mirror images. In a standard (achiral) analytical environment, they have identical physical and chemical properties. They would not be separated by a standard LC column and would produce a single, identical set of signals in a standard NMR spectrum. This is contradicted by both the LC and NMR observations. This is incorrect.
*   **D) The crude compound exists as a mixture of diastereoisomers:** Diastereoisomers are stereoisomers that are *not* mirror images of each other. They have different physical and chemical properties.
    *   Because they have different physical properties, they can be separated by standard LC, explaining the two LC peaks.
    *   Because they are chemically distinct, their corresponding protons are in different environments and will have different chemical shifts in the NMR spectrum, explaining the two NMR peaks.
    *   As isomers, they have the same mass, which is consistent with the MS data.
    *   This explanation fits all experimental observations perfectly. In peptide chemistry, this is a common phenomenon, often arising from either epimerization at one of several chiral centers or, more commonly, from restricted rotation around an amide bond (especially one involving proline), which creates stable *cis* and *trans* rotational isomers (rotamers). These rotamers are a specific type of diastereomer.

**Conclusion**

The combination of two separable peaks in LC with identical mass, along with two distinct signals for a single proton in NMR, is the classic signature of a mixture of diastereoisomers. This is the only option that is consistent with all the provided analytical data.

<<<D>>>
"""

# Run the check
result = check_chemistry_answer(llm_answer_text)
print(result)