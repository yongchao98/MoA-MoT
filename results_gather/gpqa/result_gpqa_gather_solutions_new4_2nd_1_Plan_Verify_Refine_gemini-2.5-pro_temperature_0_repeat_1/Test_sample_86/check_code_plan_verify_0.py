import re

def check_correctness():
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The logic is as follows:
    1.  **Deduce the Structure:** Based on the question's constraints (di-substituted aromatic ring, 8 total carbons, carbonyl group, aromatic-halogen bond), the only logical structure is a haloacetophenone (X-C6H4-COCH3).
    2.  **Predict the Spectrum:** For a haloacetophenone, the simplest and most common isomer is the para-substituted one. This isomer would exhibit a specific 1H NMR pattern:
        - An acetyl group (-COCH3) gives a singlet (s) for 3 protons (3H), typically around 2.1-2.7 ppm.
        - A para-substituted aromatic ring (-C6H4-) gives two doublets (d), each for 2 protons (2H), in the aromatic region (typically > 6.5 ppm).
    3.  **Evaluate Options:** The code will programmatically check each option against these predicted features.
    4.  **Compare:** The code's conclusion is compared with the LLM's provided answer.
    """

    # The final answer provided by the LLM in the prompt
    llm_answer = 'C'

    # The options from the question
    options = {
        'A': "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        'B': "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)",
        'C': "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        'D': "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)"
    }

    def parse_nmr_data(data_string):
        """Parses the NMR data string into a list of dictionaries."""
        signals = []
        # Regex to find patterns like: 7.8 (2H, d)
        pattern = re.compile(r"(\d+\.\d+)\s*\((\d+)H,\s*([a-z]+)\)")
        matches = pattern.findall(data_string)
        for match in matches:
            signals.append({
                'shift': float(match[0]),
                'integration': int(match[1]),
                'multiplicity': match[2]
            })
        return signals

    correct_option = None
    error_messages = {}

    for label, data_str in options.items():
        signals = parse_nmr_data(data_str)
        
        # --- Define expected properties for para-haloacetophenone ---
        # Structure: X-C6H4-COCH3. Total protons = 4 (aromatic) + 3 (methyl) = 7.
        
        # Check 1: Total proton count
        total_protons = sum(s['integration'] for s in signals)
        if total_protons != 7:
            error_messages[label] = f"Incorrect total proton count. Expected 7, but got {total_protons}."
            continue

        # Check 2: Presence of an aldehyde signal (which would be incorrect)
        has_aldehyde = any(s['shift'] > 9.0 for s in signals)
        if has_aldehyde:
            error_messages[label] = "Incorrect functional group. Signal > 9.0 ppm indicates an aldehyde, but the structure is a ketone."
            continue
            
        # Check 3: Aromatic protons
        aromatic_signals = [s for s in signals if s['shift'] > 6.5]
        aromatic_proton_count = sum(s['integration'] for s in aromatic_signals)
        if aromatic_proton_count != 4:
            error_messages[label] = f"Incorrect number of aromatic protons. Expected 4, but got {aromatic_proton_count}."
            continue
        
        # Check 4: Aromatic signal pattern (expecting two 2H doublets for para-substitution)
        if len(aromatic_signals) != 2 or not all(s['integration'] == 2 and s['multiplicity'] == 'd' for s in aromatic_signals):
            error_messages[label] = "Incorrect aromatic splitting pattern. Expected two 2H doublets for a para-substituted ring."
            continue

        # Check 5: Acetyl group protons
        aliphatic_signals = [s for s in signals if s['shift'] < 3.0]
        aliphatic_proton_count = sum(s['integration'] for s in aliphatic_signals)
        if aliphatic_proton_count != 3:
            error_messages[label] = f"Incorrect number of aliphatic protons. Expected 3 for an acetyl group, but got {aliphatic_proton_count}."
            continue
            
        # Check 6: Acetyl group signal pattern (expecting one 3H singlet)
        if len(aliphatic_signals) != 1 or not all(s['integration'] == 3 and s['multiplicity'] == 's' for s in aliphatic_signals):
            error_messages[label] = "Incorrect aliphatic signal pattern. Expected one 3H singlet for an acetyl group."
            continue
        
        # If all checks pass, this is the correct option
        correct_option = label
        break # Stop checking once the correct one is found

    if llm_answer == correct_option:
        return "Correct"
    else:
        if correct_option is None:
             return f"The LLM's answer '{llm_answer}' is incorrect. Furthermore, the checker could not find any option that satisfies all chemical constraints."
        
        llm_reason = error_messages.get(llm_answer, f"The LLM's choice '{llm_answer}' does not match the expected spectrum.")
        
        return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option}.\n"
                f"Reasoning:\n"
                f"The molecular structure is deduced to be para-haloacetophenone. This requires a specific 1H NMR spectrum:\n"
                f"1. Four aromatic protons (typically > 6.5 ppm), appearing as two 2H doublets.\n"
                f"2. Three aliphatic protons from an acetyl group (-COCH3), appearing as one 3H singlet (typically 2.1-2.7 ppm).\n"
                f"The LLM's choice '{llm_answer}' is incorrect because: {llm_reason}\n"
                f"Option '{correct_option}' is correct because it perfectly matches all these criteria.")

# Run the check
result = check_correctness()
print(result)