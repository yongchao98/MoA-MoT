import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by verifying the 1H NMR data against
    the constraints derived from the chemical problem description.
    """

    # The final answer provided by the LLM analysis
    llm_answer = "B"

    # The options from the question
    options = {
        "A": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "B": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "C": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        "D": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)"
    }

    def parse_nmr_data(data_string):
        """Parses an NMR data string into a list of peak dictionaries."""
        peaks = []
        pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)")
        matches = pattern.findall(data_string)
        for match in matches:
            peaks.append({
                "shift": float(match[0]),
                "integration": int(match[1]),
                "multiplicity": match[2]
            })
        return peaks

    def check_spectrum(nmr_data):
        """
        Checks if a given NMR spectrum matches the constraints for para-haloacetophenone.
        Returns "Pass" or a failure reason string.
        """
        # Constraint 1: Total proton integration must be 7 (4 aromatic + 3 methyl).
        total_protons = sum(p['integration'] for p in nmr_data)
        if total_protons != 7:
            return f"Total proton integration is {total_protons}, but should be 7."

        # Constraint 2: Must not contain an aldehyde peak (shift > 9 ppm).
        if any(p['shift'] > 9.0 for p in nmr_data):
            return "Contains a peak in the aldehyde region (>9 ppm), but the compound is a ketone."

        # Separate peaks into aromatic and aliphatic regions
        aromatic_peaks = [p for p in nmr_data if p['shift'] >= 6.5]
        aliphatic_peaks = [p for p in nmr_data if p['shift'] < 6.5]

        # Constraint 3: Must have 4 protons in the aromatic region.
        if sum(p['integration'] for p in aromatic_peaks) != 4:
            return "Does not have exactly 4 protons in the aromatic region (shift >= 6.5 ppm)."

        # Constraint 4: Aromatic protons must be two 2H doublets (para-substitution pattern).
        aromatic_doublets_2h = [p for p in aromatic_peaks if p['integration'] == 2 and p['multiplicity'] == 'd']
        if len(aromatic_doublets_2h) != 2:
            return "Aromatic region does not show the characteristic pattern of two 2H doublets."

        # Constraint 5: Must have a 3H singlet for the methyl group.
        methyl_singlet = [p for p in aliphatic_peaks if p['integration'] == 3 and p['multiplicity'] == 's']
        if len(methyl_singlet) != 1:
            return "Does not contain a 3H singlet for the acetyl methyl group."
            
        # Constraint 6: The 3H singlet must be in the typical range for a methyl ketone.
        if not (2.0 <= methyl_singlet[0]['shift'] <= 2.7):
            return f"The 3H singlet at {methyl_singlet[0]['shift']} ppm is outside the typical range (2.0-2.7 ppm) for a methyl ketone."

        return "Pass"

    # Get the data for the LLM's chosen answer
    chosen_spectrum_str = options.get(llm_answer)
    if not chosen_spectrum_str:
        return f"The provided answer '{llm_answer}' is not a valid option key."

    # Parse and check the chosen spectrum
    parsed_spectrum = parse_nmr_data(chosen_spectrum_str)
    result = check_spectrum(parsed_spectrum)

    if result == "Pass":
        return "Correct"
    else:
        return f"Incorrect. The answer '{llm_answer}' is wrong because it fails a key constraint: {result}"

# Execute the check
print(check_correctness())