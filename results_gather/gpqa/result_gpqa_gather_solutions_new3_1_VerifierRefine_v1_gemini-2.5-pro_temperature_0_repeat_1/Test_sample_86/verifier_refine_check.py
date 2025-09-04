import re

def check_correctness_of_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the given NMR spectroscopy problem.

    The function deduces the expected NMR spectrum from the question's constraints
    and checks if the chosen answer's data matches these expectations.

    Constraints from the question:
    1.  Di-substituted 6-membered aromatic ring -> C6H4 core, 4 aromatic protons.
    2.  Total 8 carbons -> Substituents have 2 carbons.
    3.  Carbonyl (C=O) and aromatic-halogen (Ar-X) -> Substituents are -X and -COCH3.
    4.  The molecule is a haloacetophenone.

    Expected 1H NMR features for the most likely isomer (para-haloacetophenone):
    -   Aromatic region: Two doublets, each integrating to 2H. Total 4H. (ppm ~7-8.5)
    -   Aliphatic region: One singlet, integrating to 3H. (ppm ~2.0-2.7)
    -   No aldehyde proton (ppm ~9-10).
    """

    # Extract the final answer choice (e.g., 'A', 'B', 'C', 'D') from the text
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Failure: Could not find a final answer in the format <<<X>>> in the provided text."
    
    chosen_option = match.group(1)

    # Define the NMR data for each option from the original question
    options = {
        'A': [{'ppm': 7.8, 'H': 2, 'split': 'd'}, {'ppm': 7.6, 'H': 2, 'split': 'd'}, {'ppm': 2.3, 'H': 3, 'split': 's'}],
        'B': [{'ppm': 4.8, 'H': 2, 'split': 'd'}, {'ppm': 4.6, 'H': 2, 'split': 'd'}, {'ppm': 1.3, 'H': 3, 'split': 's'}],
        'C': [{'ppm': 9.9, 'H': 1, 'split': 's'}, {'ppm': 7.8, 'H': 2, 'split': 'd'}, {'ppm': 7.6, 'H': 2, 'split': 'd'}, {'ppm': 3.7, 'H': 2, 'split': 's'}],
        'D': [{'ppm': 6.9, 'H': 1, 'split': 's'}, {'ppm': 4.8, 'H': 2, 'split': 'd'}, {'ppm': 4.6, 'H': 2, 'split': 'd'}, {'ppm': 1.3, 'H': 2, 'split': 's'}]
    }

    spectrum_to_check = options.get(chosen_option)
    if not spectrum_to_check:
        return f"Failure: The chosen option '{chosen_option}' is not a valid option."

    # --- Define and apply checks based on chemical principles ---

    # Check 1: Presence and integration of aromatic protons
    aromatic_protons = [p for p in spectrum_to_check if p['ppm'] >= 6.5]
    if not aromatic_protons:
        return f"Incorrect. The chosen answer '{chosen_option}' is wrong because its spectrum lacks signals in the aromatic region (ppm > 6.5), which is required for an aromatic compound."
    
    total_aromatic_H = sum(p['H'] for p in aromatic_protons)
    if total_aromatic_H != 4:
        return f"Incorrect. The chosen answer '{chosen_option}' is wrong because the total integration of its aromatic protons is {total_aromatic_H}H, but a di-substituted benzene ring must have 4H."

    # Check 2: Presence of an aldehyde signal (which should not be there)
    aldehyde_signals = [p for p in spectrum_to_check if p['ppm'] > 9.0]
    if aldehyde_signals:
        return f"Incorrect. The chosen answer '{chosen_option}' is wrong because its spectrum shows a signal above 9.0 ppm, which is characteristic of an aldehyde. The deduced structure is a ketone."

    # Check 3: Presence of the acetyl group signal (-COCH3)
    acetyl_signals = [p for p in spectrum_to_check if p['split'] == 's' and p['H'] == 3 and 2.0 <= p['ppm'] <= 2.7]
    if not acetyl_signals:
        return f"Incorrect. The chosen answer '{chosen_option}' is wrong because its spectrum does not show a singlet for 3 protons in the typical range (2.0-2.7 ppm) for an acetyl group (-COCH3)."

    # Check 4: Presence of the para-substitution pattern in the aromatic region
    aromatic_doublets_2H = [p for p in aromatic_protons if p['split'] == 'd' and p['H'] == 2]
    if len(aromatic_doublets_2H) != 2:
         return f"Incorrect. The chosen answer '{chosen_option}' is wrong because its aromatic region does not show the characteristic pattern of two doublets, each integrating to 2H, which is expected for a para-substituted ring."

    # If all checks pass, the answer is correct.
    return "Correct"

# The user's final answer to be checked is passed to the function.
final_answer_from_user = """
Here is a step-by-step analysis to determine the correct answer.

**1. Deconstruct the Information and Deduce the Molecular Structure**

*   **Core Structure:** The compound has a "di-substituted 6-membered aromatic ring," which means it is a benzene ring (C₆H₄) with two substituents.
*   **Carbon Count:** The molecule has a total of 8 carbon atoms. Since the benzene ring accounts for 6 carbons, the two substituents together must contain 2 carbon atoms (8 - 6 = 2).
*   **Functional Groups (from FTIR):**
    *   There is a **carbonyl group (C=O)**.
    *   There is an **aromatic-halogen bond**, meaning a halogen atom (F, Cl, Br, or I) is directly attached to the benzene ring.
*   **Assembling the Structure:**
    *   One substituent is a halogen atom (-X). This substituent contains 0 carbon atoms.
    *   The other substituent must therefore contain both of the remaining 2 carbon atoms and the carbonyl group. The most common and stable chemical group that fits this description is an **acetyl group (-COCH₃)**.
    *   Therefore, the compound is a **haloacetophenone** (X-C₆H₄-COCH₃). This structure perfectly matches all the given information.

**2. Predict the ¹H NMR Spectrum**

Now, we need to predict the ¹H NMR spectrum for a haloacetophenone and determine the most likely substitution pattern (ortho, meta, or para) by examining the options.

*   **Acetyl Protons (-COCH₃):** The three protons of the methyl group are adjacent to a carbonyl carbon, which has no protons. Therefore, they will not be split and will appear as a **singlet (s)**, integrating to **3H**. Their typical chemical shift is in the range of 2.0-2.7 ppm.
*   **Aromatic Protons (on the C₆H₄ ring):** There are four protons on the ring. Their splitting pattern is key:
    *   **Para (1,4) substitution:** This symmetrical arrangement results in two sets of chemically equivalent protons. This produces a very characteristic and simple pattern: **two doublets (d)**, each integrating to **2H**.
    *   **Ortho (1,2) or Meta (1,3) substitution:** These less symmetrical arrangements would result in more complex splitting patterns, typically with four distinct signals (multiplets).

Given that multiple-choice questions in spectroscopy often use symmetrical examples, the para-isomer is the most likely candidate.

**3. Evaluate the ¹H NMR Options**

Let's compare our prediction (two doublets for 4H in the aromatic region, one singlet for 3H in the aliphatic region) with the given options.

*   **A) 1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)**
    *   `7.8 (2H, d)` and `7.6 (2H, d)`: This is the classic pattern for a para-substituted benzene ring. The signals are in the correct downfield region for an aromatic ring with electron-withdrawing groups.
    *   `2.3 (3H, s)`: This is a singlet integrating to 3 protons at a chemical shift that is perfect for an acetyl methyl group.
    *   **Conclusion:** This option is an excellent match for para-haloacetophenone.

*   **B) 1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)**
    *   The signals at 4.8 and 4.6 ppm are not in the typical aromatic region (which is ~6.5-8.5 ppm). This is incorrect.

*   **C) 1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)**
    *   The signal at 9.9 ppm is characteristic of an aldehyde proton (-CHO), not a ketone. Our deduced structure is a ketone. This is incorrect.

*   **D) 1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)**
    *   This spectrum shows only one proton in the aromatic region (`6.9 (1H, s)`), which is inconsistent with a di-substituted benzene ring that must have four aromatic protons. This is incorrect.

**Final Decision**

The only option that correctly describes the ¹H NMR spectrum of the deduced molecule, para-haloacetophenone, is A. It accurately shows the singlet for the acetyl group and the pair of doublets characteristic of the para-substituted aromatic ring.

<<<A>>>
"""

# Run the check
result = check_correctness_of_answer(final_answer_from_user)
print(result)