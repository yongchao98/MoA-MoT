import re

def check_answer(question, candidate_answers, final_answer_text):
    """
    Checks the correctness of the final answer for the given NMR spectroscopy question.

    The function will:
    1.  Deduce the correct molecular structure from the question's constraints.
    2.  Predict the expected 1H NMR spectrum for that structure.
    3.  Parse the provided options (A, B, C, D) into a structured format.
    4.  Evaluate each option against the predicted spectrum.
    5.  Compare the correct option with the provided final answer.
    """

    # Step 1 & 2: Deduce structure and predict spectrum based on question constraints
    # - "di-substituted 6-membered aromatic ring": C6H4 core.
    # - "8 carbon atoms in total": 6 in ring, so 2 in substituents.
    # - "carbonyl group": C=O present.
    # - "aromatic-halogen bond": Halogen (-X) attached to the ring.
    # Deduction: One substituent is -X (0 carbons). The other must have 2 carbons and a C=O.
    # This uniquely identifies an acetyl group (-COCH3).
    # The structure is a haloacetophenone (X-C6H4-COCH3).
    #
    # Prediction for 1H NMR (assuming the most symmetrical para-isomer, which is typical for such problems):
    # - Acetyl group (-COCH3): 3 protons, next to a C=O (no H), so it's a singlet (s). Shift ~2.1-2.7 ppm.
    # - Aromatic ring (-C6H4-): 4 protons. Para-substitution gives two sets of 2 equivalent protons.
    #   This results in a characteristic pattern of two doublets (d), each for 2H. Shift ~6.5-8.5 ppm.

    # Helper function to parse the NMR data strings
    def parse_nmr_data(data_string):
        """Parses a string like '7.8 (2H, d)' into a dictionary."""
        signals = []
        # Use regex to find all occurrences of the pattern: shift (integration, splitting)
        pattern = r"(\d+\.\d+)\s*\((\d+)H,\s*([a-z]+)\)"
        matches = re.findall(pattern, data_string)
        for match in matches:
            signals.append({
                'shift': float(match[0]),
                'protons': int(match[1]),
                'splitting': match[2]
            })
        return signals

    # Step 3: Parse the options from the question text
    # The options are consistent in the question prompt itself.
    options_text = {
        'A': "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        'B': "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        'C': "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        'D': "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)"
    }
    
    parsed_options = {key: parse_nmr_data(value) for key, value in options_text.items()}

    # Step 4: Evaluate each option
    correct_option = None
    for key, signals in parsed_options.items():
        # Reset checks for each option
        has_methyl_ketone_signal = False
        aromatic_protons = 0
        aromatic_doublets_2h = 0
        has_aldehyde_proton = False
        total_protons = sum(s['protons'] for s in signals)

        for signal in signals:
            # Check for the acetyl methyl group signal
            if 2.1 <= signal['shift'] <= 2.7 and signal['protons'] == 3 and signal['splitting'] == 's':
                has_methyl_ketone_signal = True
            
            # Check for aromatic signals
            if signal['shift'] >= 6.5:
                # Check for aldehyde proton (common distractor)
                if signal['shift'] > 9.0:
                    has_aldehyde_proton = True
                else: # It's an aromatic proton
                    aromatic_protons += signal['protons']
                    if signal['protons'] == 2 and signal['splitting'] == 'd':
                        aromatic_doublets_2h += 1
        
        # A spectrum is correct if it meets all criteria for para-haloacetophenone
        if (has_methyl_ketone_signal and
            not has_aldehyde_proton and
            aromatic_protons == 4 and
            aromatic_doublets_2h == 2 and
            len(signals) == 3): # Expecting exactly 3 signals
            correct_option = key
            break

    # Step 5: Compare with the final answer
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Could not find a final answer in the standard format '<<<X>>>' in the provided text."
    
    final_answer_key = match.group(1)

    if final_answer_key == correct_option:
        return "Correct"
    else:
        # Provide a detailed reason for the error
        if correct_option is None:
             return f"The provided answer '{final_answer_key}' is incorrect. In fact, none of the options perfectly match the expected spectrum for a haloacetophenone."

        # Analyze why the chosen answer is wrong
        chosen_option_data = parsed_options.get(final_answer_key, [])
        reason = ""
        if not chosen_option_data:
            reason = f"The chosen answer '{final_answer_key}' is not a valid option."
        else:
            # Example checks for common errors
            if any(s['shift'] > 9.0 for s in chosen_option_data):
                reason = "it contains a signal > 9.0 ppm, which is characteristic of an aldehyde, not the required ketone."
            elif sum(s['protons'] for s in chosen_option_data if s['shift'] > 6.5) != 4:
                reason = "it does not account for the 4 protons on the di-substituted aromatic ring."
            elif not any(2.1 <= s['shift'] <= 2.7 and s['protons'] == 3 and s['splitting'] == 's' for s in chosen_option_data):
                reason = "it lacks the characteristic 3H singlet for the acetyl group (-COCH3)."
            else:
                reason = "it does not fit the expected pattern for a haloacetophenone."

        return f"Incorrect. The provided answer '{final_answer_key}' is wrong because {reason} The correct answer is '{correct_option}', which shows the characteristic pattern for para-haloacetophenone: two 2H doublets in the aromatic region (7.8 and 7.6 ppm) and one 3H singlet for the acetyl group (2.3 ppm)."

# Extract the question and final answer from the provided text
question_text = """
You have prepared a di-substituted 6-membered aromatic ring compound containing 8 carbon atoms in total. Peaks in the FTIR spectrum of this compound indicate the presence of a carbonyl group as well as an aromatic-halogen bond. Identify the correct 1H NMR data for this compound.

A) 1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)
B) 1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)
C) 1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)
D) 1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)
"""

final_answer_text_from_llm = """
Here is a step-by-step analysis to determine the final answer, taking into account the provided information and the reasoning from the candidate answers.

### Step 1: Deduce the Molecular Structure

First, the chemical structure of the compound must be determined from the information given in the question.

*   **"di-substituted 6-membered aromatic ring"**: The core of the molecule is a benzene ring with two substituents. The ring itself accounts for 6 carbon atoms and has the formula -C₆H₄-.
*   **"8 carbon atoms in total"**: Since the benzene ring has 6 carbons, the two substituents together must contain the remaining 2 carbon atoms.
*   **"FTIR spectrum...indicate the presence of a carbonyl group"**: One of the substituents contains a C=O group.
*   **"FTIR spectrum...indicate the presence of an aromatic-halogen bond"**: One substituent is a halogen atom (F, Cl, Br, or I) attached directly to the aromatic ring.

To assemble the structure, we combine these facts:
1.  One substituent is a halogen (-X). It contains 0 carbon atoms.
2.  The other substituent must contain the remaining 2 carbons and the carbonyl group. The most common and stable chemical group that fits this description is an **acetyl group (-COCH₃)**.

Therefore, the compound is a **haloacetophenone**, with the general structure X-C₆H₄-COCH₃. This structure satisfies all the initial conditions.

### Step 2: Predict the ¹H NMR Spectrum

Next, we predict the ¹H NMR spectrum for a haloacetophenone. The appearance of the spectrum depends on the relative positions of the two substituents (ortho, meta, or para).

*   **Acetyl Group Protons (-COCH₃)**: The three protons of the methyl group are adjacent to a carbonyl carbon, which has no protons. Therefore, they will not be split and will appear as a **singlet (s)**, integrating to **3H**. The typical chemical shift for these protons is in the range of δ 2.1-2.7 ppm.
*   **Aromatic Protons (-C₆H₄-)**: There are four protons on the ring.
    *   If the substitution is **para (1,4)**, the molecule has a plane of symmetry. This results in a simple and highly characteristic pattern: **two doublets (d)**, each integrating to **2H**.
    *   If the substitution is ortho (1,2) or meta (1,3), the molecule is asymmetrical, which would lead to a more complex spectrum with up to four distinct signals for the four non-equivalent protons.

Given that the options provided are all simple patterns, it is overwhelmingly likely that the compound is the **para isomer**. Our prediction is therefore a spectrum containing a 3H singlet and two 2H doublets.

### Step 3: Evaluate the Given Options

Now, we must carefully compare our prediction to the specific options listed in the original question prompt.

*   **A) 1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)**
    *   This is incorrect. It shows only one aromatic proton (at 6.9 ppm), but a di-substituted ring must have four.

*   **B) 1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)**
    *   This is incorrect. The signal at 9.9 ppm is characteristic of an aldehyde proton (-CHO), not the ketone (-COCH₃) that was deduced from the structural information.

*   **C) 1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)**
    *   This is incorrect. It completely lacks signals in the typical aromatic region (δ > 6.5 ppm), which is impossible for an aromatic compound.

*   **D) 1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)**
    *   This is a **perfect match** for our prediction of para-haloacetophenone.
        *   `2.3 (3H, s)`: This is the singlet for the three protons of the acetyl methyl group.
        *   `7.8 (2H, d)` and `7.6 (2H, d)`: This is the classic pattern of two doublets for the four protons on a para-substituted aromatic ring.

### Careful Point: Analysis of Candidate Answers

The provided candidate answers almost unanimously perform the correct chemical analysis, identifying the structure as para-haloacetophenone and the correct NMR data as being two aromatic doublets and one methyl singlet. However, many of them incorrectly map this correct data to the letters A, B, or C due to an apparent shuffling of the options. The sound chemical reasoning across all answers points to the data set presented in option D of the original question.

### Conclusion

The molecular structure is haloacetophenone, and the NMR data corresponds to the para-isomer. Based on a direct comparison with the options provided in the question, only option D accurately represents the ¹H NMR spectrum for this compound.

<<<D>>>
"""

# Run the check
result = check_answer(question_text, [], final_answer_text_from_llm)
print(result)