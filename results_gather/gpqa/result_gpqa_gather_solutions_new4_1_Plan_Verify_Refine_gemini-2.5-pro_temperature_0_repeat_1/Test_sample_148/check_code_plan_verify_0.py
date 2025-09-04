import re

def check_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer based on the problem's constraints.

    The function analyzes the experimental data provided in the question:
    1.  ¹H NMR: Two peaks for one proton -> two distinct species, distinguishable by NMR.
    2.  LC-MS: Two peaks -> two species, separable by LC.
    3.  LC-MS: Same mass for both peaks -> the two species are isomers.

    It then evaluates the chosen answer against these facts.
    """
    # Extract the final answer letter (A, B, C, or D) from the LLM's response.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<A>>>, <<<B>>>, <<<C>>>, or <<<D>>>."
    
    final_answer = match.group(1)

    # Define the key facts derived from the experimental data.
    # Fact 1: The two species have the same mass (they are isomers).
    # This is from the LC-MS data showing both peaks have the same mass spectrum.
    species_are_isomers = True

    # Fact 2: The two species are distinguishable by standard NMR.
    # This is from the ¹H NMR data showing two distinct peaks for the same proton.
    species_are_distinguishable_by_nmr = True

    # Fact 3: The two species are separable by standard LC.
    # This is from the LC data showing two clearly defined peaks.
    species_are_separable_by_lc = True

    # Define the properties of each possible explanation (the options).
    explanations = {
        'A': {  # 'Double coupling'
            "description": "'Double coupling' product",
            "is_isomer": False,  # Would have a different (higher) mass.
            "is_distinguishable_by_nmr": True,
            "is_separable_by_lc": True,
        },
        'B': {  # Precursor contamination
            "description": "Precursor contamination",
            "is_isomer": False,  # Would have a different mass.
            "is_distinguishable_by_nmr": True,
            "is_separable_by_lc": True,
        },
        'C': {  # Enantiomers
            "description": "Mixture of enantiomers",
            "is_isomer": True,
            "is_distinguishable_by_nmr": False,  # Not distinguishable by standard NMR.
            "is_separable_by_lc": False,  # Not separable by standard LC.
        },
        'D': {  # Diastereomers
            "description": "Mixture of diastereoisomers",
            "is_isomer": True,
            "is_distinguishable_by_nmr": True,
            "is_separable_by_lc": True,
        }
    }

    # Get the properties of the chosen answer.
    chosen_explanation = explanations.get(final_answer)

    # Check if the chosen explanation's properties match the experimental facts.
    if chosen_explanation["is_isomer"] != species_are_isomers:
        return (f"Incorrect. The answer '{final_answer}' ({chosen_explanation['description']}) is wrong. "
                f"The LC-MS data shows both species have the same mass, meaning they must be isomers. "
                f"A {chosen_explanation['description']} would have a different molecular weight.")

    if chosen_explanation["is_separable_by_lc"] != species_are_separable_by_lc:
        return (f"Incorrect. The answer '{final_answer}' ({chosen_explanation['description']}) is wrong. "
                f"The LC data shows two distinct peaks, meaning the species are separable by standard chromatography. "
                f"Enantiomers are not separable by standard LC.")

    if chosen_explanation["is_distinguishable_by_nmr"] != species_are_distinguishable_by_nmr:
        return (f"Incorrect. The answer '{final_answer}' ({chosen_explanation['description']}) is wrong. "
                f"The ¹H NMR data shows two distinct signals for the same proton, meaning the species are distinguishable by standard NMR. "
                f"Enantiomers are not distinguishable by standard NMR.")

    # If all properties match the facts, the answer is correct.
    return "Correct"

# The provided answer from the LLM
llm_answer = """
Here is a step-by-step analysis to determine the final answer.

**Plan:**
*   Systematically analyze the key pieces of information from the ¹H NMR and LC-MS data.
*   Evaluate each of the four provided options against the combined experimental evidence.
*   Synthesize the findings to identify the single most plausible explanation.

**Execution:**

**1. Analysis of Experimental Data:**
*   **¹H NMR Data:** The observation of two distinct peaks for a single alpha-proton, with equal integrals summing to 1H, indicates that the proton exists in two different chemical environments. The two molecular species containing this proton are present in a 1:1 ratio.
*   **LC-MS Data:** The liquid chromatography (LC) separates the sample into two distinct peaks of equal intensity, confirming the presence of two physically separable species in a 1:1 ratio. The mass spectrometry (MS) data is critical: both peaks yield the same mass spectrum, which matches the expected product. This proves that the two species are **isomers** (same chemical formula, same mass) and not impurities with different molecular weights.

**2. Evaluation of Options:**

*   **A) 'Double coupling' has occurred during an amide-bond forming reaction:** This side reaction would create a product with a different, higher molecular weight. This is directly contradicted by the MS data showing both peaks have the *expected* mass.
*   **B) The compound is contaminated with a precursor:** A precursor would be a different molecule with a different molecular weight. This is also contradicted by the MS data.
*   **C) The crude compound exists as a mixture of enantiomers:** Enantiomers are non-superimposable mirror images. In a standard (achiral) analytical environment, they have identical physical properties. They would produce a single set of signals in the NMR and a single peak in the LC. This contradicts the observation of two peaks in both experiments.
*   **D) The crude compound exists as a mixture of diastereoisomers:** Diastereoisomers are stereoisomers that are *not* mirror images. They have different physical and chemical properties.
    *   **Consistency with NMR:** Because they are chemically distinct, their corresponding protons are in different environments and will have different chemical shifts, explaining the two NMR peaks.
    *   **Consistency with LC-MS:** Because they have different physical properties, they can be separated by standard chromatography, explaining the two LC peaks. As isomers, they have the same mass, which is consistent with the MS data. This is a classic signature of a diastereomeric mixture, often caused in peptides by epimerization or by *cis-trans* isomerism around an amide bond (atropisomerism).

**3. Synthesis:**
The only explanation that is consistent with all pieces of experimental data is that the compound exists as a mixture of diastereoisomers. The candidate answers, despite some confusion in their final letter choices, overwhelmingly support this conclusion in their reasoning.

<<<D>>>
"""

# Run the check and print the result.
result = check_correctness(llm_answer)
print(result)