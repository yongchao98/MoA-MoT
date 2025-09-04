import re

def check_answer():
    """
    Checks the correctness of the given answer for the 1H NMR spectroscopy question.
    """
    question_constraints = {
        "ring_type": "di-substituted 6-membered aromatic ring",
        "total_carbons": 8,
        "has_carbonyl": True,
        "has_aromatic_halogen": True
    }

    # The final answer provided by the LLM
    final_answer = "C"

    options = {
        "A": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "B": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "C": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "D": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)"
    }

    # Step 1: Deduce the structure and expected NMR features from the question constraints.
    # - Aromatic ring (6C) + 2 substituents (2C) = 8C total.
    # - Substituents: Halogen (-X) and a 2C group with a carbonyl. This must be an acetyl group (-COCH3).
    # - The structure is a haloacetophenone (X-C6H4-COCH3).
    # - Total protons = 4 on the ring + 3 on the methyl group = 7 protons.
    # - Expected NMR for the most symmetrical isomer (para-haloacetophenone), which matches the simple patterns in the options:
    #   1. Aromatic protons: 4H total, appearing as two doublets (2H each) in the 7.0-8.0 ppm range.
    #   2. Acetyl methyl protons: 3H, appearing as a singlet in the 2.1-2.6 ppm range.
    
    selected_option_data = options.get(final_answer)
    if not selected_option_data:
        return f"Invalid answer option '{final_answer}'. Please choose from A, B, C, D."

    # Use regex to parse the NMR data string
    # Pattern matches: ppm_value (integrationH, splitting)
    pattern = re.compile(r"(\d+\.\d+)\s*\((\d+)H,\s*([a-z]+)\)")
    signals = pattern.findall(selected_option_data)

    if not signals:
        return f"Could not parse the NMR data for option {final_answer}."

    # Convert parsed strings to numerical types
    parsed_signals = []
    for ppm, integration, splitting in signals:
        parsed_signals.append({
            "ppm": float(ppm),
            "H": int(integration),
            "split": splitting
        })

    # Step 2: Check the selected answer against the expected features.
    
    # Check 1: Total proton count
    total_protons = sum(s['H'] for s in parsed_signals)
    if total_protons != 7:
        return f"Incorrect. The total proton count for a haloacetophenone (C8H7XO) should be 7, but the answer has {total_protons} protons."

    # Check 2: Aromatic protons
    aromatic_signals = [s for s in parsed_signals if s['ppm'] >= 6.5 and s['split'] == 'd' and s['H'] == 2]
    if len(aromatic_signals) != 2:
        return f"Incorrect. The spectrum should show two doublets for 2H each in the aromatic region (characteristic of para-substitution). The selected answer does not meet this criterion."

    # Check 3: Acetyl methyl protons
    acetyl_signal = [s for s in parsed_signals if 2.0 <= s['ppm'] <= 2.7 and s['split'] == 's' and s['H'] == 3]
    if len(acetyl_signal) != 1:
        return f"Incorrect. The spectrum should show a singlet for 3H in the 2.0-2.7 ppm range for the acetyl methyl group. This signal is missing or incorrect in the selected answer."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)