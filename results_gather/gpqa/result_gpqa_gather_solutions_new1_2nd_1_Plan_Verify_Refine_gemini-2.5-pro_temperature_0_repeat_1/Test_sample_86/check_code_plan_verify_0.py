import re

def check_correctness():
    """
    This function checks the correctness of the provided LLM answer by analyzing the chemical constraints of the question.

    The question describes a di-substituted 8-carbon aromatic ring with a carbonyl and an aromatic-halogen bond.
    This uniquely points to a haloacetophenone (C8H7XO).
    The simple NMR patterns in the options suggest a para-substituted isomer.

    The expected 1H NMR for para-haloacetophenone is:
    1. Aromatic protons: 4H total, appearing as two 2H doublets in the 6.5-8.5 ppm range.
    2. Methyl protons: A 3H singlet for the acetyl group (-COCH3) in the 2.0-2.7 ppm range.
    3. Total protons: 7H.
    4. No aldehyde proton (which would be > 9 ppm).
    """

    # The options as presented in the original question context
    options = {
        "A": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "B": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "C": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "D": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)"
    }

    # The final answer provided by the LLM
    llm_answer = "C"

    def parse_nmr_data(data_string):
        """Parses a string of NMR data into a list of dictionaries."""
        peaks = []
        pattern = re.compile(r"(\d+\.?\d+)\s*\((\d+)H,\s*([a-z]+)\)")
        peak_strings = data_string.split(',')
        for peak_str in peak_strings:
            match = pattern.search(peak_str.strip())
            if match:
                peaks.append({
                    "ppm": float(match.group(1)),
                    "integration": int(match.group(2)),
                    "splitting": match.group(3).strip()
                })
        return peaks

    def check_spectrum(peaks):
        """Checks if a spectrum matches the expected properties. Returns (is_match, reason)."""
        # Constraint 1: Total proton count must be 7 (for C8H7XO).
        total_protons = sum(p['integration'] for p in peaks)
        if total_protons != 7:
            return False, f"Fails total proton count constraint. Expected 7H for C8H7XO, but found {total_protons}H."

        # Constraint 2: Must not be an aldehyde.
        if any(p['ppm'] >= 9.0 for p in peaks):
            return False, "Fails functional group constraint. Contains an aldehyde peak (ppm >= 9.0), but the compound is a ketone."

        # Constraint 3: Must have aromatic protons.
        aromatic_peaks = [p for p in peaks if 6.5 <= p['ppm'] < 9.0]
        if not aromatic_peaks:
            return False, "Fails aromatic constraint. No signals found in the aromatic region (6.5-9.0 ppm)."

        # Constraint 4: Aromatic protons must match para-substitution pattern (4H total, two 2H doublets).
        aromatic_proton_count = sum(p['integration'] for p in aromatic_peaks)
        if aromatic_proton_count != 4:
            return False, f"Fails aromatic proton count. Expected 4H, but found {aromatic_proton_count}H."
        
        aromatic_doublets = [p for p in aromatic_peaks if p['splitting'] == 'd' and p['integration'] == 2]
        if len(aromatic_doublets) != 2 or len(aromatic_peaks) != 2:
            return False, "Fails aromatic splitting pattern. Expected two 2H doublets for para-substitution."

        # Constraint 5: Must have a methyl ketone signal (3H singlet).
        methyl_ketone_peaks = [p for p in peaks if 2.0 <= p['ppm'] <= 2.7 and p['integration'] == 3 and p['splitting'] == 's']
        if len(methyl_ketone_peaks) != 1:
            return False, "Fails methyl ketone constraint. Expected a 3H singlet between 2.0-2.7 ppm."

        return True, "Matches all constraints for para-haloacetophenone."

    # Check the LLM's chosen answer
    selected_option_data = parse_nmr_data(options[llm_answer])
    is_correct, reason = check_spectrum(selected_option_data)

    if not is_correct:
        return f"The answer '{llm_answer}' is incorrect. Reason: {reason}"

    # Verify that all other options are incorrect to ensure the answer is unique
    for key, value in options.items():
        if key == llm_answer:
            continue
        other_option_data = parse_nmr_data(value)
        is_other_correct, _ = check_spectrum(other_option_data)
        if is_other_correct:
            return f"The answer '{llm_answer}' is plausible, but option '{key}' also satisfies all constraints. The question may be ambiguous or the check is not specific enough."

    return "Correct"

# Run the check and print the result
print(check_correctness())