import re

def check_nmr_correctness():
    """
    This function checks the correctness of the provided answer by analyzing the chemical constraints
    of the question and applying them to the selected 1H NMR data.

    The function will:
    1.  Deduce the expected molecular structure and its corresponding 1H NMR features.
    2.  Define the options given in the question and the final answer provided.
    3.  Check if the NMR data of the chosen answer satisfies all the deduced chemical constraints.
    4.  Return "Correct" if it matches, or a reason for the mismatch if it doesn't.
    """

    # Step 1: Deduce chemical structure and expected 1H NMR features from the question.
    # - Di-substituted 6-membered aromatic ring: A benzene ring with two groups (C6H4).
    # - 8 carbon atoms total: Ring (6C) + Substituents (2C).
    # - FTIR shows carbonyl (C=O) and aromatic-halogen (Ar-X).
    # - Deduction: One substituent is a halogen (-X, 0 carbons). The other must be an acetyl group (-COCH3, 2 carbons, 1 C=O).
    # - Structure: Haloacetophenone (X-C6H4-COCH3).
    # - The simple splitting patterns in the options (doublets) strongly suggest para-substitution.
    # - Expected 1H NMR for para-haloacetophenone (C8H7XO):
    #   - Total Protons: 7H.
    #   - Aromatic Protons: 4H total, in the aromatic region (>6.5 ppm), split into two 2H doublets.
    #   - Acetyl Protons: 3H, as a singlet, in the typical methyl ketone region (~2.0-2.7 ppm).
    #   - No aldehyde proton (which would appear > 9.0 ppm).

    # Step 2: Define the options from the question and the LLM's final answer.
    options = {
        "A": "6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "B": "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "C": "4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "D": "9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)"
    }
    llm_final_answer = "B"

    # Step 3: Get the NMR data for the chosen answer.
    chosen_nmr_data_string = options.get(llm_final_answer)
    if not chosen_nmr_data_string:
        return f"Invalid answer choice '{llm_final_answer}'. Not one of the options A, B, C, D."

    # Step 4: Parse the NMR data string and check against constraints.
    signals = []
    # Regex to capture ppm (float), integration (int), and multiplicity (string)
    pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)")
    matches = pattern.findall(chosen_nmr_data_string)
    
    if not matches:
        return f"Could not parse the NMR data string: '{chosen_nmr_data_string}'"
        
    for match in matches:
        signals.append({
            "ppm": float(match[0]),
            "integration": int(match[1]),
            "multiplicity": match[2]
        })

    # Constraint Check 1: Disqualifying features (e.g., aldehyde proton).
    if any(s['ppm'] > 9.0 for s in signals):
        aldehyde_ppm = next(s['ppm'] for s in signals if s['ppm'] > 9.0)
        return f"Incorrect: The answer contains a signal at {aldehyde_ppm} ppm, which is characteristic of an aldehyde. The deduced structure is a ketone."

    # Constraint Check 2: Total proton count.
    total_protons = sum(s['integration'] for s in signals)
    if total_protons != 7:
        return f"Incorrect total proton count. Expected 7 for C8H7XO, but the data sums to {total_protons}H."

    # Constraint Check 3: Acetyl group signal (-COCH3).
    acetyl_signals = [s for s in signals if s['multiplicity'] == 's' and s['integration'] == 3]
    if len(acetyl_signals) != 1:
        return "Incorrect: The data does not contain the required 3H singlet for the acetyl group."
    acetyl_ppm = acetyl_signals[0]['ppm']
    if not (2.0 <= acetyl_ppm <= 2.7):
        return f"Incorrect: The 3H singlet at {acetyl_ppm} ppm is outside the expected chemical shift range (2.0-2.7 ppm) for an acetyl group."

    # Constraint Check 4: Aromatic proton signals.
    aromatic_signals = [s for s in signals if s['ppm'] > 6.5]
    aromatic_proton_count = sum(s['integration'] for s in aromatic_signals)
    if aromatic_proton_count != 4:
        return f"Incorrect number of aromatic protons. Expected 4, but found {aromatic_proton_count} in the region > 6.5 ppm."

    # Constraint Check 5: Aromatic splitting pattern (para-substitution).
    aromatic_doublets = [s for s in aromatic_signals if s['multiplicity'] == 'd' and s['integration'] == 2]
    if len(aromatic_doublets) != 2:
        return "Incorrect aromatic splitting pattern. Expected two 2H doublets for para-substitution, but this pattern was not found."

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_nmr_correctness()
print(result)