import re

def check_answer(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The function will:
    1.  Deduce the expected molecular properties from the question.
    2.  Deduce the expected 1H NMR spectral features for the molecule.
    3.  Parse the NMR data from the options.
    4.  Check if the chosen answer's data matches the expected features.
    5.  Check if the other options are correctly ruled out.
    """

    # Step 1 & 2: Deduce expected structure and NMR features from the question
    # - Di-substituted 6-membered aromatic ring -> Benzene ring with 4H
    # - 8 carbon atoms total -> 6 in ring, 2 in substituents
    # - Carbonyl group (C=O) and aromatic-halogen bond -> Substituents are -X and -COCH3
    # - The molecule is a haloacetophenone.
    # - Expected 1H NMR for the most likely isomer (para-haloacetophenone):
    #   - 4 aromatic protons (on the C6H4 ring)
    #   - 3 protons from the acetyl methyl group (-COCH3)
    #   - Total protons = 7
    #   - Acetyl group: 3H singlet, chemical shift ~2.1-2.7 ppm.
    #   - Aromatic protons: 4H total, chemical shift > 6.5 ppm. For para-substitution, this is typically two 2H doublets.
    #   - No aldehyde proton (shift ~9-10 ppm).

    # Define expected features for the correct answer
    expected_features = {
        "total_aromatic_protons": 4,
        "has_acetyl_group": True,
        "has_aldehyde_group": False,
        "is_aromatic": True,
        "is_para_substituted": True
    }

    # Parse the options from the question text (assuming a standard format)
    # Note: The provided LLM answers sometimes have different labels (A, B, C, D) for the same data.
    # We will define the data sets and check them against the logic.
    data_sets = {
        "set1": "9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "set2": "4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "set3": "6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "set4": "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)"
    }
    
    # The question has options A, B, C, D. We map them to the sets.
    # This mapping is consistent with the final provided answer block.
    options_map = {
        "A": data_sets["set1"],
        "B": data_sets["set2"],
        "C": data_sets["set3"],
        "D": data_sets["set4"]
    }

    def parse_nmr_data(data_string):
        """Parses a string of NMR data into a list of peak dictionaries."""
        peaks = []
        # Regex to find: chemical shift, integration, multiplicity
        pattern = re.compile(r"(\d+\.?\d+)\s*\((\d+)H,\s*([a-z])\)")
        matches = pattern.findall(data_string)
        for match in matches:
            peaks.append({
                "shift": float(match[0]),
                "integration": int(match[1]),
                "multiplicity": match[2]
            })
        return peaks

    def analyze_peaks(peaks):
        """Analyzes a list of peaks and returns a dictionary of features."""
        analysis = {
            "total_aromatic_protons": 0,
            "has_acetyl_group": False,
            "has_aldehyde_group": False,
            "is_aromatic": False,
            "is_para_substituted": False,
            "aromatic_peaks": []
        }

        for peak in peaks:
            # Check for aromatic protons
            if peak["shift"] > 6.5:
                analysis["is_aromatic"] = True
                analysis["total_aromatic_protons"] += peak["integration"]
                analysis["aromatic_peaks"].append(peak)

            # Check for acetyl group (-COCH3)
            if 2.1 <= peak["shift"] <= 2.7 and peak["integration"] == 3 and peak["multiplicity"] == 's':
                analysis["has_acetyl_group"] = True

            # Check for aldehyde group (-CHO)
            if 9.0 <= peak["shift"] <= 10.0 and peak["integration"] == 1 and peak["multiplicity"] == 's':
                analysis["has_aldehyde_group"] = True
        
        # Check for para-substitution pattern (two 2H doublets in aromatic region)
        aromatic_peaks = analysis["aromatic_peaks"]
        if len(aromatic_peaks) == 2:
            peak1, peak2 = aromatic_peaks[0], aromatic_peaks[1]
            if (peak1["integration"] == 2 and peak1["multiplicity"] == 'd' and
                peak2["integration"] == 2 and peak2["multiplicity"] == 'd'):
                analysis["is_para_substituted"] = True

        return analysis

    # Find the final answer chosen by the LLM
    final_answer_match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not final_answer_match:
        return "Could not find a final answer in the format <<<A>>>, <<<B>>>, etc."
    
    chosen_option = final_answer_match.group(1)
    
    # Get the data for the chosen option
    chosen_data_string = options_map.get(chosen_option)
    if not chosen_data_string:
        return f"The chosen option '{chosen_option}' is invalid."

    # Analyze the chosen data
    chosen_analysis = analyze_peaks(parse_nmr_data(chosen_data_string))

    # Check against expected features
    if not chosen_analysis["is_aromatic"]:
        return f"Incorrect. The chosen answer {chosen_option} is not an aromatic compound as it lacks signals above 6.5 ppm, but the question specifies an aromatic ring."
    
    if chosen_analysis["total_aromatic_protons"] != expected_features["total_aromatic_protons"]:
        return f"Incorrect. The chosen answer {chosen_option} has {chosen_analysis['total_aromatic_protons']} aromatic protons, but a di-substituted benzene ring must have 4."

    if chosen_analysis["has_aldehyde_group"]:
        return f"Incorrect. The chosen answer {chosen_option} shows a signal for an aldehyde proton (~9.9 ppm), but the deduced structure is a ketone (haloacetophenone)."

    if not chosen_analysis["has_acetyl_group"]:
        return f"Incorrect. The chosen answer {chosen_option} does not show a 3H singlet around 2.1-2.7 ppm, which is characteristic of the required acetyl group (-COCH3)."

    if not chosen_analysis["is_para_substituted"]:
         return f"Incorrect. The chosen answer {chosen_option} does not show the characteristic two 2H doublets for a para-substituted aromatic ring."

    # If all checks pass for the chosen answer, it is correct.
    # We can also double-check that the other options fail.
    correct_option = None
    for option, data_string in options_map.items():
        analysis = analyze_peaks(parse_nmr_data(data_string))
        is_correct = all(analysis[key] == expected_features[key] for key in expected_features)
        if is_correct:
            correct_option = option
            break
    
    if chosen_option == correct_option:
        return "Correct"
    else:
        return f"Incorrect. The LLM chose {chosen_option}, but the only option that satisfies all constraints is {correct_option}."


# You would pass the question and the LLM's full response text to this function.
# For this example, we'll use the provided text.
question_text = """
**Task**: Decide the final answer based on the following answers from other agents.

**Question**:
You have prepared a di-substituted 6-membered aromatic ring compound containing 8 carbon atoms in total. Peaks in the FTIR spectrum of this compound indicate the presence of a carbonyl group as well as an aromatic-halogen bond. Identify the correct 1H NMR data for this compound.

A) 1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)
B) 1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)
C) 1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)
D) 1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)
"""

llm_response = """
Here is a step-by-step analysis of the question to determine the correct answer.

**Step 1: Deduce the Molecular Structure**

1.  **Core Structure:** The compound is described as a "di-substituted 6-membered aromatic ring." This means the base of the molecule is a benzene ring with two substituents attached. The ring itself accounts for 6 carbon atoms and has the formula -C₆H₄-.
2.  **Carbon Count:** The molecule has a total of 8 carbon atoms. Since the benzene ring has 6 carbons, the two substituents together must contain the remaining 2 carbon atoms (8 - 6 = 2).
3.  **Functional Groups (from FTIR):**
    *   The presence of a "carbonyl group" (C=O) is indicated.
    *   The presence of an "aromatic-halogen bond" is indicated, meaning one substituent is a halogen atom (-X, where X = F, Cl, Br, or I) attached directly to the ring.
4.  **Assembling the Structure:**
    *   One substituent is a halogen (-X), which contains no carbon atoms.
    *   Therefore, the other substituent must contain both of the remaining 2 carbons and the carbonyl group.
    *   The most common and stable chemical group that fits this description is an **acetyl group (-COCH₃)**.
    *   Thus, the compound is a **haloacetophenone**, with the general structure X-C₆H₄-COCH₃. This structure satisfies all the given conditions.

**Step 2: Predict the ¹H NMR Spectrum**

Now, we predict the ¹H NMR spectrum for a haloacetophenone. The spectrum's appearance depends on the relative positions of the two substituents (ortho, meta, or para). The simple splitting patterns presented in the options strongly suggest the most symmetrical isomer, the **para (1,4) isomer**.

*   **Acetyl Group Protons (-COCH₃):** The three protons of the methyl group are adjacent to a carbonyl carbon, which has no protons. Therefore, they will not be split and will appear as a **singlet (s)**, integrating to **3H**. Their chemical shift is typically in the range of δ 2.1–2.7 ppm.
*   **Aromatic Protons (-C₆H₄-):** In a para-substituted ring, there is a plane of symmetry. This results in two sets of chemically equivalent protons. This leads to a very characteristic and simple pattern: **two doublets (d)**, each integrating to **2H**. These signals will appear in the downfield aromatic region (typically δ 7.0–8.5 ppm) due to the electron-withdrawing effects of the halogen and the acetyl group.

**Step 3: Evaluate the Given Options**

We will now compare our prediction for *para*-haloacetophenone with the provided choices. We are looking for a spectrum with two 2H doublets in the aromatic region and one 3H singlet around 2.3 ppm.

*   **A) 1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)**
    *   This is incorrect. The signal at δ 9.9 ppm is characteristic of an **aldehyde** proton (-CHO), not the ketone (-COCH₃) that our structure contains.

*   **B) 1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)**
    *   This is incorrect. It lacks any signals in the typical aromatic region (usually δ > 6.5 ppm), which is required for an aromatic compound.

*   **C) 1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)**
    *   This is incorrect. It shows only one proton in the aromatic region (`6.9 (1H, s)`), which is inconsistent with a di-substituted benzene ring that must have four aromatic protons.

*   **D) 1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)**
    *   This is a **perfect match** for our prediction.
        *   `7.8 (2H, d)` and `7.6 (2H, d)`: This is the classic pattern for the four protons on a para-substituted aromatic ring.
        *   `2.3 (3H, s)`: This is the singlet for the three protons of the acetyl methyl group, with the correct integration and a typical chemical shift.

**Conclusion:**
The molecular structure is a para-haloacetophenone. The ¹H NMR data in option D is the only one that accurately represents this structure.

<<<D>>>
"""

result = check_answer(question_text, llm_response)
print(result)