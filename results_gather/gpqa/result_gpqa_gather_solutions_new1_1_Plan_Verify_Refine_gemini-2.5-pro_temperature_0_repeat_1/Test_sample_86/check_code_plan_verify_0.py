import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by analyzing the chemical constraints.
    """
    # The final answer provided by the LLM
    llm_answer = 'D'

    # The options from the question
    options = {
        'A': "6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        'B': "9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        'C': "4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        'D': "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)"
    }

    def parse_nmr_data(data_string):
        """Parses an NMR data string into a list of dictionaries."""
        peaks = []
        pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)")
        matches = pattern.findall(data_string)
        for match in matches:
            peaks.append({
                'ppm': float(match[0]),
                'H': int(match[1]),
                'split': match[2].strip()
            })
        return peaks

    def check_spectrum(spectrum_data):
        """
        Checks if a spectrum matches the expected data for a para-haloacetophenone.
        The structure X-C6H4-COCH3 has 7 protons.
        """
        # Constraint 1: Total number of protons must be 7.
        total_protons = sum(p['H'] for p in spectrum_data)
        if total_protons != 7:
            return False, f"Incorrect total proton count. Expected 7, found {total_protons}."

        # Constraint 2: Must have a signal for the acetyl methyl group.
        # This is a 3H singlet, typically between 2.0 and 3.0 ppm.
        methyl_signal = [p for p in spectrum_data if p['H'] == 3 and p['split'] == 's' and 2.0 <= p['ppm'] <= 3.0]
        if len(methyl_signal) != 1:
            return False, "Does not contain the characteristic 3H singlet for an acetyl group in the expected region."

        # Constraint 3: Must have 4 protons in the aromatic region (6.5-8.5 ppm).
        aromatic_signals = [p for p in spectrum_data if 6.5 <= p['ppm'] <= 8.5]
        total_aromatic_protons = sum(p['H'] for p in aromatic_signals)
        if total_aromatic_protons != 4:
            return False, f"Incorrect number of aromatic protons. Expected 4, found {total_aromatic_protons}."

        # Constraint 4: The aromatic signals should match a para-substitution pattern (two 2H doublets).
        para_pattern_signals = [p for p in aromatic_signals if p['H'] == 2 and p['split'] == 'd']
        if len(para_pattern_signals) != 2:
            return False, "Aromatic signals do not match the pattern for para-substitution (two 2H doublets)."

        return True, "Spectrum is consistent with para-haloacetophenone."

    # Find which option is chemically correct
    correct_option_key = None
    for key, data in options.items():
        spectrum = parse_nmr_data(data)
        is_valid, reason = check_spectrum(spectrum)
        if is_valid:
            correct_option_key = key
            break

    if correct_option_key is None:
        return "Checker failed: No option satisfies all the chemical constraints."

    if llm_answer == correct_option_key:
        return "Correct"
    else:
        is_llm_answer_valid, reason_llm = check_spectrum(parse_nmr_data(options[llm_answer]))
        return (f"Incorrect. The provided answer is {llm_answer}, but the chemically correct answer is {correct_option_key}.\n"
                f"Reasoning for why {llm_answer} is wrong: {reason_llm}")

# Execute the check
result = check_correctness()
print(result)