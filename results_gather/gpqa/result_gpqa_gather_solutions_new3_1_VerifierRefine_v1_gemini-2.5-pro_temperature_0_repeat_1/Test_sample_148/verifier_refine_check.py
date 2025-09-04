import re

def check_peptide_isomer_answer(llm_response: str):
    """
    Checks the correctness of the answer to the peptide analysis question.

    The function simulates the logical deduction process based on the provided
    NMR and LC-MS data.

    Args:
        llm_response: The full response from the LLM, containing the answer
                      in the format <<<ANSWER>>>.

    Returns:
        A string indicating "Correct" or a detailed reason for the incorrectness.
    """
    # --- Step 1: Define the facts from the experimental data ---
    # Fact 1: LC-MS shows two peaks with the SAME mass, consistent with the expected molecule.
    # This means the two species are ISOMERS of the expected molecule.
    facts = {
        "are_isomers": True,
        # Fact 2: 1H NMR shows two distinct peaks for one proton.
        # This means the isomers are distinguishable by standard NMR.
        "distinguishable_by_nmr": True,
        # Fact 3: LC shows two clearly defined peaks.
        # This means the isomers are separable by standard LC.
        "separable_by_lc": True
    }

    # --- Step 2: Define the properties of each option from the question ---
    # A) 'Double coupling' -> Higher mass, not an isomer.
    # B) Enantiomers -> Isomers, but NOT distinguishable/separable in achiral conditions.
    # C) Diastereoisomers -> Isomers, and ARE distinguishable/separable.
    # D) Precursor -> Lower mass, not an isomer.
    options_properties = {
        'A': {"name": "'Double coupling' product", "is_isomer": False, "distinguishable_by_nmr": True, "separable_by_lc": True},
        'B': {"name": "Enantiomers", "is_isomer": True, "distinguishable_by_nmr": False, "separable_by_lc": False},
        'C': {"name": "Diastereoisomers", "is_isomer": True, "distinguishable_by_nmr": True, "separable_by_lc": True},
        'D': {"name": "Precursor", "is_isomer": False, "distinguishable_by_nmr": True, "separable_by_lc": True}
    }

    # --- Step 3: Find the logically correct answer ---
    logically_correct_option = None
    for option, properties in options_properties.items():
        if (properties["is_isomer"] == facts["are_isomers"] and
            properties["distinguishable_by_nmr"] == facts["distinguishable_by_nmr"] and
            properties["separable_by_lc"] == facts["separable_by_lc"]):
            logically_correct_option = option
            break

    # --- Step 4: Extract and validate the LLM's answer ---
    match = re.search(r'<<<([A-D])>>>', llm_response)
    if not match:
        return "The provided answer is not in the correct format '<<<X>>>' or is not a valid option."
    
    llm_answer = match.group(1)

    # --- Step 5: Compare and generate the final result ---
    if llm_answer == logically_correct_option:
        return "Correct"
    else:
        chosen_properties = options_properties[llm_answer]
        reason = f"The answer '{llm_answer}' ({chosen_properties['name']}) is incorrect.\n"
        
        # Check the mass/isomer constraint
        if chosen_properties['is_isomer'] != facts['are_isomers']:
            reason += "Constraint failed: The LC-MS data shows both peaks have the same mass as the expected molecule, meaning they must be isomers. "
            reason += f"A {chosen_properties['name']} is not an isomer of the target molecule."
            return reason

        # Check the distinguishability/separability constraint
        if chosen_properties['distinguishable_by_nmr'] != facts['distinguishable_by_nmr']:
            reason += "Constraint failed: The NMR shows two distinct peaks, but enantiomers are indistinguishable in a standard (achiral) NMR experiment."
            return reason
            
        if chosen_properties['separable_by_lc'] != facts['separable_by_lc']:
            reason += "Constraint failed: The LC shows two separate peaks, but enantiomers are not separable on a standard (achiral) LC column."
            return reason
            
        return f"An unknown error occurred in reasoning for answer {llm_answer}."


# The final response from the LLM to be checked
llm_final_response = """
Here is a step-by-step analysis of the question to determine the final answer.

**Step 1: Analyze the Mass Spectrometry (MS) Data**
The LC-MS analysis shows two peaks, but crucially, both peaks have the **same mass spectrum**, and this mass is consistent with the **expected molecule**.
*   **Implication:** The two species present in the sample are isomers. They have the same molecular formula and the same mass.
*   **Elimination:** This observation immediately rules out options that would result in a different mass.
    *   **D) The compound is contaminated with a precursor:** A precursor is an incomplete molecule and would have a lower mass. This is incorrect.
    *   **A) 'Double coupling' has occurred during an amide-bond forming reaction:** This side reaction would involve adding an extra amino acid, resulting in a product with a higher mass. This is incorrect.

**Step 2: Analyze the ¹H NMR and Liquid Chromatography (LC) Data**
The ¹H NMR spectrum shows **two distinct peaks** for a single proton, and the LC analysis shows **two clearly defined peaks**.
*   **Implication:** The two isomers are physically and spectroscopically distinct in a standard (achiral) environment. They can be separated by chromatography, and they produce different signals in an NMR spectrum.

**Step 3: Evaluate the Remaining Isomer Options (Enantiomers vs. Diastereoisomers)**
We have established that the sample is a mixture of isomers. The remaining options are B) enantiomers and C) diastereoisomers.
*   **B) The crude compound exists as a mixture of enantiomers:** Enantiomers are non-superimposable mirror images. In an achiral environment (like a standard NMR solvent or a standard LC column), they have identical physical and chemical properties. Therefore, a mixture of enantiomers would produce only a single set of peaks in the NMR and would co-elute as a single peak in the LC. This contradicts the experimental observations.
*   **C) The crude compound exists as a mixture of diastereoisomers:** Diastereoisomers are stereoisomers that are *not* mirror images of each other. They have different physical and chemical properties. Because they are distinct compounds, they are expected to have different NMR spectra (explaining the two peaks for the alpha-proton) and can be separated by standard chromatography (explaining the two LC peaks). This explanation is perfectly consistent with all the data.

**Conclusion**
The evidence collectively points to a single conclusion: the sample is a mixture of diastereoisomers. They are isomers (same mass), but they are not mirror images, which allows them to be distinguished and separated by standard NMR and LC techniques.

<<<C>>>
"""

# Run the check
result = check_peptide_isomer_answer(llm_final_response)
print(result)