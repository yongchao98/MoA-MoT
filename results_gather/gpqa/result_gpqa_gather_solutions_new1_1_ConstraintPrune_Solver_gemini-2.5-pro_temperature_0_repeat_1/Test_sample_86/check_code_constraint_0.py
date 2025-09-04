import re

def check_correctness_of_nmr_answer():
    """
    Checks the correctness of the LLM's answer for the given chemistry problem.
    
    The problem describes a di-substituted 8-carbon aromatic ring with a carbonyl
    and an aromatic-halogen bond. This points to a haloacetophenone structure.
    The code checks if the NMR data for the chosen answer matches this structure.
    """

    # --- Problem Constraints & Expected NMR Features ---
    # Deduced structure: Haloacetophenone (C8H7XO)
    # Expected features for the most likely para-isomer:
    expected_total_protons = 7
    expected_aromatic_protons = 4
    expected_acetyl_protons = 3
    aromatic_ppm_range = (6.5, 8.5)
    acetyl_ppm_range = (2.0, 2.7)
    aldehyde_ppm_range = (9.5, 10.5)

    # --- Data from the Question ---
    options = {
        'A': "6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        'B': "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        'C': "9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        'D': "4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer = "B"

    def parse_nmr_data(data_string):
        """Parses a string of NMR data into a list of dictionaries."""
        signals = []
        parts = re.findall(r'([\d\.]+) \((\d+)H, ([a-z]+)\)', data_string)
        for part in parts:
            signals.append({
                'ppm': float(part[0]),
                'integration': int(part[1]),
                'splitting': part[2]
            })
        return signals

    # --- Verification Logic ---
    if llm_answer not in options:
        return f"Invalid answer option '{llm_answer}'. Please choose from A, B, C, D."

    data_string = options[llm_answer]
    parsed_data = parse_nmr_data(data_string)

    # Check 1: Total number of protons
    total_protons = sum(s['integration'] for s in parsed_data)
    if total_protons != expected_total_protons:
        return f"Incorrect. The total number of protons in option {llm_answer} is {total_protons}, but it should be {expected_total_protons} for a haloacetophenone (C8H7XO)."

    # Check 2: Presence and number of aromatic protons
    aromatic_signals = [s for s in parsed_data if aromatic_ppm_range[0] <= s['ppm'] <= aromatic_ppm_range[1]]
    total_aromatic_protons = sum(s['integration'] for s in aromatic_signals)
    if total_aromatic_protons != expected_aromatic_protons:
        return f"Incorrect. The total number of aromatic protons (in the {aromatic_ppm_range} ppm range) is {total_aromatic_protons}, but it should be {expected_aromatic_protons} for a di-substituted benzene ring."

    # Check 3: Presence of an acetyl group signal (-COCH3)
    acetyl_signal = [s for s in parsed_data if 
                     acetyl_ppm_range[0] <= s['ppm'] <= acetyl_ppm_range[1] and
                     s['integration'] == expected_acetyl_protons and
                     s['splitting'] == 's']
    if not acetyl_signal:
        return f"Incorrect. The spectrum lacks the characteristic signal for an acetyl group (a 3H singlet in the {acetyl_ppm_range} ppm range)."

    # Check 4: Absence of aldehyde protons
    aldehyde_signal = [s for s in parsed_data if aldehyde_ppm_range[0] <= s['ppm'] <= aldehyde_ppm_range[1]]
    if aldehyde_signal:
        return f"Incorrect. The spectrum shows a signal at {aldehyde_signal[0]['ppm']} ppm, which is characteristic of an aldehyde, not the expected ketone."

    # Check 5: Specific pattern for para-substitution
    aromatic_doublets = [s for s in aromatic_signals if s['splitting'] == 'd' and s['integration'] == 2]
    if len(aromatic_doublets) != 2 or len(aromatic_signals) != 2:
        return f"Incorrect. The aromatic region does not show the classic pattern of two 2H doublets expected for a para-substituted ring, which is the most plausible structure."

    return "Correct"

# Execute the check and print the result
result = check_correctness_of_nmr_answer()
print(result)