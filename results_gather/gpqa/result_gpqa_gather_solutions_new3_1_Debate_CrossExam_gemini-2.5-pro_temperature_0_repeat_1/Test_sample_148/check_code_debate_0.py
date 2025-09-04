import re

def check_peptide_analysis_answer(llm_response: str) -> str:
    """
    Checks the correctness of the answer to the peptide analysis question.

    The function models the question's constraints as a set of logical rules
    and checks if the provided answer satisfies all of them.
    """
    try:
        # Extract the final answer choice (e.g., 'D') from the response text.
        match = re.search(r'<<<([A-D])>>>', llm_response)
        if not match:
            return "Incorrect: The answer format is wrong or missing. It should be <<<X>>> where X is A, B, C, or D."
        final_answer = match.group(1)
    except Exception:
        return "Incorrect: Could not parse the final answer from the response."

    # 1. Define the experimental facts derived from the question's text.
    # - LC-MS shows two peaks with the SAME mass as the expected molecule -> The components are isomers.
    # - ¹H NMR shows two distinct peaks for a single proton -> The components are distinguishable by standard NMR.
    # - LC shows two separable peaks -> The components are separable by standard (achiral) chromatography.
    experimental_facts = {
        "are_isomers": True,
        "distinguishable_by_nmr": True,
        "separable_by_lc": True
    }

    # 2. Define the theoretical properties of each possible explanation.
    # The options are mapped as follows:
    # A: Precursor
    # B: Double coupling
    # C: Enantiomers
    # D: Diastereoisomers
    properties = {
        "A": {  # Precursor
            "are_isomers": False,  # A precursor has a different (lower) mass.
            "distinguishable_by_nmr": True,
            "separable_by_lc": True
        },
        "B": {  # 'Double coupling' product
            "are_isomers": False,  # This product would have a different (higher) mass.
            "distinguishable_by_nmr": True,
            "separable_by_lc": True
        },
        "C": {  # Enantiomers
            "are_isomers": True,
            "distinguishable_by_nmr": False,  # Indistinguishable in a standard (achiral) NMR.
            "separable_by_lc": False  # Not separable by a standard (achiral) LC column.
        },
        "D": {  # Diastereoisomers
            "are_isomers": True,
            "distinguishable_by_nmr": True,
            "separable_by_lc": True
        }
    }

    # 3. Compare the properties of the chosen answer with the experimental facts.
    chosen_properties = properties.get(final_answer)
    
    mismatches = []
    for key, fact_value in experimental_facts.items():
        if chosen_properties[key] != fact_value:
            mismatches.append(key)

    # 4. Return the verdict.
    if not mismatches:
        return "Correct"
    else:
        reasons = []
        for mismatch in mismatches:
            if mismatch == "are_isomers":
                reasons.append("it is inconsistent with the MS data, which shows both components have the same mass (i.e., they are isomers). A precursor or double-coupled product would have a different mass.")
            elif mismatch == "distinguishable_by_nmr":
                reasons.append("it is inconsistent with the NMR data, which shows two distinct peaks for one proton, meaning the components are distinguishable in a standard NMR experiment. Enantiomers are not.")
            elif mismatch == "separable_by_lc":
                reasons.append("it is inconsistent with the LC data, which shows two separate peaks, meaning the components are separable by standard chromatography. Enantiomers are not.")
        
        return f"Incorrect: The answer '{final_answer}' is wrong because " + " and ".join(reasons)

# The provided LLM response to be checked.
llm_answer_text = """
### Self-Debate

**Proposer:** The consensus among all detailed candidate answers is clear and points to a single explanation. The analytical data—specifically two separable LC peaks with identical mass, and two distinct NMR signals for a single proton—is the classic signature of a mixture of diastereoisomers. The identical mass from the MS data definitively rules out precursors (Option A) and 'double coupling' products (Option B), as both would have different molecular weights. The separability by standard LC and distinguishability by standard NMR rules out enantiomers (Option C), which are identical in achiral environments. This leaves only diastereoisomers (Option D) as the viable explanation.

**Skeptic:** The logic is sound, and the process of elimination is convincing. However, the problem mentions the LC-MS was run at an "elevated temperature." Some answers suggest this is to improve separation, while others note that high temperatures can cause rapidly interconverting isomers (like some rotamers) to coalesce into a single peak. Since the peaks remained separate, does this imply we are dealing with a specific type of stable diastereomer, like one formed from epimerization, rather than cis-trans rotamers?

**Proposer:** That's a fair point to consider, but it doesn't change the final answer. Both epimers and cis-trans rotamers are types of diastereomers. The energy barrier for cis-trans isomerization around a peptide bond, especially one involving proline, can be quite high, often requiring significant thermal energy to cause rapid interconversion on the NMR or LC timescale. The fact that the peaks remain separate at an elevated temperature simply confirms that under these specific analytical conditions, the two species are stable and not rapidly equilibrating. It reinforces that we are observing two distinct chemical entities. Whether they are epimers or rotamers, they both fall under the umbrella of diastereomers. Therefore, the conclusion that the sample is a mixture of diastereoisomers remains the most accurate and encompassing explanation for all the data presented.

**Skeptic:** I agree. The core evidence is that the two species are isomers (same mass) but have different physical properties (separable by LC) and different spectroscopic properties (distinct NMR signals). This is the textbook definition of how diastereomers behave in an analytical setting. The other options are fundamentally inconsistent with the data. The debate converges.

---
### Final Analysis

Based on a step-by-step analysis of the provided experimental data, the correct conclusion is reached through a process of elimination:

1.  **Analyze the LC-MS Data:** The crude compound separates into two distinct peaks in the liquid chromatogram (LC). The mass spectrometer (MS) shows that both peaks have the *same mass*, which is consistent with the expected product.
    *   **Conclusion:** The sample contains a mixture of two **isomers** (compounds with the same chemical formula but different structures).
    *   This observation immediately **rules out Option A (precursor)** and **Option B (double coupling)**, as both would result in products with a different molecular weight.

2.  **Analyze the ¹H NMR and LC Data Together:** The ¹H NMR spectrum shows two distinct peaks for a single proton, and the LC trace shows two separable peaks.
    *   **Conclusion:** The two isomers have different chemical and physical properties, as they can be distinguished by NMR and separated by chromatography under standard (achiral) conditions.
    *   **Enantiomers** (Option C) are mirror-image stereoisomers that have identical properties in an achiral environment. They would produce a single set of peaks in a standard NMR and co-elute as a single peak in a standard LC. This contradicts the observations.
    *   **Diastereoisomers** (Option D) are stereoisomers that are *not* mirror images. They have different physical and chemical properties. Consequently, they give distinct NMR spectra and can be separated by standard chromatography. This perfectly matches all the experimental evidence.

This phenomenon is common in peptide chemistry, often arising from either epimerization at a chiral center or, more frequently, from the slow cis-trans isomerization around a peptide bond (particularly with proline), which creates stable, separable diastereomeric conformers called atropisomers.

<<<D>>>
"""

# Run the check
result = check_peptide_analysis_answer(llm_answer_text)
print(result)