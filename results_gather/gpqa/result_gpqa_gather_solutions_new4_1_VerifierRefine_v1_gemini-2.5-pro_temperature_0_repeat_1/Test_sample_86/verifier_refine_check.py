import re

def parse_nmr_data(data_string):
    """Parses a string of NMR data into a list of dictionaries."""
    signals = []
    # Example format: "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)"
    pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-zA-Z]+)\)")
    matches = pattern.findall(data_string)
    for match in matches:
        signals.append({
            "shift": float(match[0]),
            "integration": int(match[1]),
            "multiplicity": match[2]
        })
    return signals

def check_correctness():
    """
    Checks the correctness of the answer based on the question's constraints.
    """
    # --- Step 1: Define constraints from the question ---
    # The structure is a di-substituted 6-membered aromatic ring with 8 total carbons,
    # a carbonyl group, and an aromatic-halogen bond.
    # This deduces the structure to be a haloacetophenone (X-C6H4-COCH3).
    # The most likely isomer for a simple spectrum is the para-isomer.

    # --- Step 2: Define the expected NMR spectrum for para-haloacetophenone ---
    expected_total_protons = 7  # 4 on the ring, 3 on the methyl group
    
    # Expected features for the acetyl methyl group
    expected_methyl = {
        "integration": 3,
        "multiplicity": "s",
        "shift_min": 2.1,
        "shift_max": 2.7
    }

    # Expected features for the aromatic protons (para-substituted)
    expected_aromatic_protons = 4
    expected_aromatic_shift_min = 6.5
    expected_aromatic_shift_max = 8.5

    # --- Step 3: Define the options and the given answer ---
    options = {
        "A": "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)",
        "B": "1H NMR: chemical reference (ppm): 6.9 (1H, s), 4.8 (2H, d), 4.6 (2H, d), 1.3 (2H, s)",
        "C": "1H NMR: chemical reference (ppm): 9.9 (1H, s), 7.8 (2H, d), 7.6 (2H, d), 3.7 (2H, s)",
        "D": "1H NMR: chemical reference (ppm): 4.8 (2H, d), 4.6 (2H, d), 1.3 (3H, s)"
    }
    
    given_answer_key = "A"
    answer_data_str = options[given_answer_key]
    
    # --- Step 4: Parse and verify the given answer ---
    parsed_signals = parse_nmr_data(answer_data_str)

    # Check 1: Total proton count
    total_integration = sum(s['integration'] for s in parsed_signals)
    if total_integration != expected_total_protons:
        return f"Incorrect. The total number of protons should be {expected_total_protons} (4 aromatic + 3 methyl), but the answer has a total integration of {total_integration}H."

    # Check 2: Presence and properties of the acetyl methyl signal
    methyl_signal_found = False
    for signal in parsed_signals:
        if (signal['integration'] == expected_methyl['integration'] and
            signal['multiplicity'] == expected_methyl['multiplicity'] and
            expected_methyl['shift_min'] <= signal['shift'] <= expected_methyl['shift_max']):
            methyl_signal_found = True
            break
    if not methyl_signal_found:
        return f"Incorrect. The spectrum is missing the characteristic signal for an acetyl methyl group (a 3H singlet between {expected_methyl['shift_min']}-{expected_methyl['shift_max']} ppm)."

    # Check 3: Presence and properties of the aromatic signals
    aromatic_protons_count = 0
    aromatic_signals = []
    for signal in parsed_signals:
        if expected_aromatic_shift_min <= signal['shift'] <= expected_aromatic_shift_max:
            aromatic_protons_count += signal['integration']
            aromatic_signals.append(signal)
            
    if aromatic_protons_count != expected_aromatic_protons:
        return f"Incorrect. A di-substituted benzene ring must have {expected_aromatic_protons} aromatic protons, but the answer only accounts for {aromatic_protons_count} in the aromatic region ({expected_aromatic_shift_min}-{expected_aromatic_shift_max} ppm)."

    # Check 4: Specific pattern for para-substitution (most likely case)
    # We expect two doublets, each for 2H.
    doublet_2h_count = 0
    for signal in aromatic_signals:
        if signal['integration'] == 2 and signal['multiplicity'] == 'd':
            doublet_2h_count += 1
    
    if doublet_2h_count != 2:
         return f"Incorrect. The pattern in the aromatic region does not match the classic 'two 2H doublets' expected for a para-substituted ring. Found {doublet_2h_count} signals matching this pattern."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)