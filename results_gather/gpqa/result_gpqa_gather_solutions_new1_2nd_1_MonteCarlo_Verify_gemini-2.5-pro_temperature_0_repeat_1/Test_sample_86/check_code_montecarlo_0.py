import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer by deducing the compound's structure
    and verifying its expected 1H NMR spectrum against the given options.
    """

    # --- 1. Define the problem and the LLM's answer ---
    question_options = {
        'A': "6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        'B': "9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        'C': "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        'D': "4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
    }
    llm_answer = 'C'

    # --- 2. Helper function to parse NMR data strings ---
    def parse_nmr_data(data_string):
        """Parses a string of NMR data into a list of signal dictionaries."""
        signals = []
        # Regex to find patterns like: 7.8 (2H, d)
        pattern = re.compile(r'(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)')
        matches = pattern.findall(data_string)
        for match in matches:
            signals.append({
                'shift': float(match[0]),
                'integration': int(match[1]),
                'multiplicity': match[2]
            })
        return signals

    parsed_options = {key: parse_nmr_data(value) for key, value in question_options.items()}

    # --- 3. Define chemical constraint checks based on the deduced structure (para-haloacetophenone) ---

    def check_total_protons(signals):
        # Structure is C8H7XO, so total protons must be 7.
        total_h = sum(s['integration'] for s in signals)
        if total_h != 7:
            return False, f"Fails proton count: Expected 7 total protons, but found {total_h}."
        return True, ""

    def check_aromatic_protons(signals):
        # A di-substituted benzene ring must have 4 aromatic protons (shift > 6.5 ppm).
        aromatic_h = sum(s['integration'] for s in signals if s['shift'] > 6.5)
        if aromatic_h != 4:
            return False, f"Fails aromatic proton count: Expected 4 protons with shift > 6.5 ppm, but found {aromatic_h}."
        return True, ""

    def check_ketone_not_aldehyde(signals):
        # Structure is a ketone, not an aldehyde. No signal should be in the aldehyde region (> 9.0 ppm).
        if any(s['shift'] > 9.0 for s in signals):
            return False, "Fails ketone check: Found a signal > 9.0 ppm, which indicates an aldehyde, not a ketone."
        return True, ""

    def check_acetyl_group(signals):
        # Structure has an acetyl group (-COCH3), which should be a 3H singlet.
        # Typical shift is ~2.0-2.7 ppm.
        if not any(s['integration'] == 3 and s['multiplicity'] == 's' for s in signals):
            return False, "Fails acetyl group check: Did not find the required 3H singlet."
        return True, ""
        
    def check_para_substitution(signals):
        # The simple pattern suggests para-substitution: two 2H doublets in the aromatic region.
        aromatic_signals = [s for s in signals if s['shift'] > 6.5]
        if len(aromatic_signals) != 2:
            return False, f"Fails para-substitution check: Expected 2 signals in the aromatic region, but found {len(aromatic_signals)}."
        
        doublet_2h_count = sum(1 for s in aromatic_signals if s['integration'] == 2 and s['multiplicity'] == 'd')
        if doublet_2h_count != 2:
            return False, "Fails para-substitution check: The aromatic signals are not two 2H doublets."
        return True, ""

    # --- 4. Evaluate the LLM's chosen answer ---
    
    # All checks that must be passed by the correct answer
    all_checks = [
        check_total_protons,
        check_aromatic_protons,
        check_ketone_not_aldehyde,
        check_acetyl_group,
        check_para_substitution
    ]

    signals_to_check = parsed_options.get(llm_answer)
    if not signals_to_check:
        return f"Invalid answer key '{llm_answer}'. It is not one of the options A, B, C, D."

    for check_func in all_checks:
        passed, reason = check_func(signals_to_check)
        if not passed:
            return f"Incorrect. The chosen answer {llm_answer} is wrong because it {reason}"

    # --- 5. Verify that no other option is also correct ---
    for option_key, signals in parsed_options.items():
        if option_key == llm_answer:
            continue
        
        is_also_correct = True
        for check_func in all_checks:
            passed, _ = check_func(signals)
            if not passed:
                is_also_correct = False
                break
        
        if is_also_correct:
            return f"Incorrect. The analysis is flawed because option {option_key} also satisfies all constraints, but the final answer was {llm_answer}."

    return "Correct"

# Execute the check and print the result
print(check_correctness())