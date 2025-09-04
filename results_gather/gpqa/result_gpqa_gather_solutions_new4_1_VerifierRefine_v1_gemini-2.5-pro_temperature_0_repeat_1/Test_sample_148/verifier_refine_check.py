import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of an LLM's answer to a chemistry problem.

    The function simulates the logical deduction process based on the
    experimental evidence provided in the question.

    Args:
        llm_answer_text: The full text of the LLM's answer, including the
                         final answer in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or a detailed reason for the error.
    """

    # 1. Define the key experimental observations from the question.
    observations = {
        "num_nmr_peaks_for_1H": 2,
        "num_lc_peaks": 2,
        "mass_of_lc_peaks": "same_and_expected",  # Can be 'same_and_expected', 'different', or 'higher'
        "analysis_environment": "achiral"  # Standard NMR/LC is achiral
    }

    # 2. Define the expected analytical results for each possible option.
    option_properties = {
        "A": {  # 'Double coupling'
            "explanation": "This side reaction creates a product with a higher molecular weight.",
            "expected_mass": "higher",
            "expected_nmr_peaks": "irrelevant",
            "expected_lc_peaks": "irrelevant"
        },
        "B": {  # Contaminated with a precursor
            "explanation": "A precursor would have a different molecular weight.",
            "expected_mass": "different",
            "expected_nmr_peaks": "irrelevant",
            "expected_lc_peaks": "irrelevant"
        },
        "C": {  # Mixture of enantiomers
            "explanation": "Enantiomers are indistinguishable in an achiral environment, giving a single peak in standard NMR and LC.",
            "expected_mass": "same_and_expected",
            "expected_nmr_peaks": 1,
            "expected_lc_peaks": 1
        },
        "D": {  # Mixture of diastereoisomers
            "explanation": "Diastereoisomers are distinct compounds with different properties, leading to separate peaks in NMR and LC, but they have the same mass.",
            "expected_mass": "same_and_expected",
            "expected_nmr_peaks": ">1",
            "expected_lc_peaks": ">1"
        }
    }

    # 3. Determine the logically correct option based on the observations.
    determined_correct_option = None
    for option, properties in option_properties.items():
        # First, check if the mass observation is consistent with the option.
        mass_consistent = (properties["expected_mass"] == observations["mass_of_lc_peaks"])
        if not mass_consistent:
            continue  # Rule out this option if mass doesn't match.

        # If mass is consistent, check peak counts for isomeric options (C and D).
        if properties["expected_mass"] == "same_and_expected":
            nmr_peaks_consistent = (observations["num_nmr_peaks_for_1H"] > 1 if properties["expected_nmr_peaks"] == ">1" else observations["num_nmr_peaks_for_1H"] == properties["expected_nmr_peaks"])
            lc_peaks_consistent = (observations["num_lc_peaks"] > 1 if properties["expected_lc_peaks"] == ">1" else observations["num_lc_peaks"] == properties["expected_lc_peaks"])

            if nmr_peaks_consistent and lc_peaks_consistent:
                determined_correct_option = option
                break
        else:
            # For non-isomeric options (A and B), mass is the only key check.
            determined_correct_option = option
            break

    # 4. Extract the final answer from the LLM's response.
    # Note: The provided answer has options A, B, C, D re-mapped in some cases.
    # We will check the final letter provided.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Could not find a final answer in the format <<<A>>>, <<<B>>>, etc. in the provided text."
    
    proposed_answer = match.group(1)

    # 5. Compare the proposed answer with the determined correct answer and provide feedback.
    # The question's options are: A) Double coupling, B) Precursor, C) Enantiomers, D) Diastereomers.
    # The provided answer uses these letters.
    if proposed_answer == determined_correct_option:
        return "Correct"
    else:
        reason = f"The proposed answer '{proposed_answer}' is incorrect.\n"
        
        # Explain why the proposed answer is wrong based on its properties.
        proposed_props = option_properties[proposed_answer]
        reason += f"\nExplanation for why '{proposed_answer}' is wrong: {proposed_props['explanation']}\n"
        
        # Check mass consistency for the proposed answer
        if proposed_props["expected_mass"] != observations["mass_of_lc_peaks"]:
            reason += f"- The MS data shows the mass is '{observations['mass_of_lc_peaks']}', but option '{proposed_answer}' implies the mass would be '{proposed_props['expected_mass']}'. This is a contradiction.\n"
        
        # Check peak consistency for the proposed answer
        if "expected_nmr_peaks" in proposed_props and proposed_props["expected_nmr_peaks"] != "irrelevant":
             if not (observations["num_nmr_peaks_for_1H"] > 1 if proposed_props["expected_nmr_peaks"] == ">1" else observations["num_nmr_peaks_for_1H"] == proposed_props["expected_nmr_peaks"]):
                 reason += f"- The NMR data shows {observations['num_nmr_peaks_for_1H']} peaks for the proton, but option '{proposed_answer}' would result in {proposed_props['expected_nmr_peaks']} peak(s). This is a contradiction.\n"
        
        if "expected_lc_peaks" in proposed_props and proposed_props["expected_lc_peaks"] != "irrelevant":
             if not (observations["num_lc_peaks"] > 1 if proposed_props["expected_lc_peaks"] == ">1" else observations["num_lc_peaks"] == proposed_props["expected_lc_peaks"]):
                 reason += f"- The LC data shows {observations['num_lc_peaks']} peaks, but option '{proposed_answer}' would result in {proposed_props['expected_lc_peaks']} peak(s). This is a contradiction.\n"

        # Explain why the correct answer is right
        correct_props = option_properties[determined_correct_option]
        reason += f"\nThe correct answer is '{determined_correct_option}' because it is the only option consistent with all observations: {correct_props['explanation']}"
        
        return reason

# The provided answer to check
llm_answer = """
Here is a step-by-step analysis of the question and the provided answers.

**Step 1: Deconstruct the Experimental Observations**

1.  **¹H NMR Data:**
    *   **Observation:** Two distinct peaks for a single alpha-proton.
    *   **Implication:** This proton exists in two different chemical environments. Since spin-spin coupling is ruled out, this means there are two distinct molecular species present in the sample.
    *   **Observation:** The peaks have roughly equal integrals that sum to 1H.
    *   **Implication:** The two molecular species are present in a roughly 1:1 ratio.

2.  **LC-MS Data:**
    *   **Observation:** Two clearly defined peaks in the liquid chromatogram (LC).
    *   **Implication:** The two molecular species have different physical properties, allowing them to be separated by chromatography. This confirms they are different compounds.
    *   **Observation:** Both LC peaks have the same mass spectrum, which is consistent with the expected molecule.
    *   **Implication:** This is a critical point. It proves that the two species are **isomers**—they have the same molecular formula and mass but a different arrangement of atoms. This rules out impurities that would have a different mass.

**Step 2: Evaluate the Options Based on the Evidence**

Let's systematically check each option against our deductions from the data.

*   **A) 'Double coupling' has occurred during an amide-bond forming reaction:** This side reaction would involve adding an extra amino acid, creating a product with a different, higher molecular weight. This is contradicted by the mass spectrometry data, which shows both species have the same, expected mass. This option is incorrect.

*   **B) The compound is contaminated with a precursor:** A precursor is a starting material or intermediate. It would have a different chemical formula and therefore a different molecular weight than the final product. This is also contradicted by the mass spectrometry data. This option is incorrect.

*   **C) The crude compound exists as a mixture of enantiomers:** Enantiomers are stereoisomers that are non-superimposable mirror images. In a standard (achiral) environment, they have identical physical and chemical properties.
    *   They would not be separated by a standard LC column (they would elute as a single peak).
    *   They would not produce different signals in a standard ¹H NMR spectrum (they would show a single set of peaks).
    *   This option is inconsistent with both the LC and NMR data, which clearly show two distinct species. This option is incorrect.

*   **D) The crude compound exists as a mixture of diastereoisomers:** Diastereoisomers are stereoisomers that are *not* mirror images of each other.
    *   **Physical Properties:** They have different physical properties. This means they can be separated by standard chromatography techniques like LC, which explains the two observed LC peaks.
    *   **Chemical Properties:** They are distinct chemical compounds. This means their corresponding protons are in different chemical environments and will have different chemical shifts in an NMR spectrum, which explains the two observed NMR signals for the alpha-proton.
    *   **Mass:** As isomers, they have the same molecular formula and thus the same mass. This is consistent with both LC peaks giving the same mass spectrum.
    *   This option is the only one that perfectly explains all the experimental observations. In peptide chemistry, this is a common phenomenon, often arising from either epimerization at a chiral center or, more commonly, from the presence of stable *cis/trans* isomers (rotamers) around a peptide bond, which are a specific type of diastereomer.

**Step 3: Final Conclusion**

The combination of two separable peaks in LC with identical mass, along with two distinct signals for a single proton in NMR, is the classic signature of a mixture of diastereoisomers.

<<<D>>>
"""

# Run the check
result = check_answer_correctness(llm_answer)
print(result)