import re

def parse_nmr_data(data_string):
    """Parses a string of NMR data into a list of signal dictionaries."""
    signals = []
    # Split the string by commas that are followed by a space and a number
    parts = re.split(r',\s*(?=\d)', data_string)
    for part in parts:
        try:
            match = re.match(r'\s*(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)', part.strip())
            if match:
                shift = float(match.group(1))
                integration = int(match.group(2))
                multiplicity = match.group(3)
                signals.append({
                    "shift": shift,
                    "integration": integration,
                    "multiplicity": multiplicity
                })
            else:
                # Handle cases where parsing might fail for an invalid format
                return None
        except (ValueError, IndexError):
            return None
    return signals

def check_answer():
    """
    Checks the correctness of the provided answer based on the question's constraints.
    """
    # 1. Define constraints from the question
    # Structure: Di-substituted 6-membered aromatic ring with 8 total carbons.
    # FTIR: Has a carbonyl (C=O) and an aromatic-halogen bond.
    # This deduces the structure to be a haloacetophenone (X-C6H4-COCH3).

    # 2. Define expected properties based on the deduced structure
    expected_total_protons = 7  # 4 on the ring, 3 on the methyl group
    expected_aromatic_protons = 4
    expected_methyl_protons = 3
    
    # The most likely isomer for a simple spectrum is para-substituted.
    # Expected signals for para-haloacetophenone:
    # - Aromatic: Two 2H doublets in the range ~7.0-8.5 ppm.
    # - Aliphatic: One 3H singlet for the methyl ketone in the range ~2.0-2.7 ppm.

    # 3. The final answer to check is D
    answer_choice = 'D'
    answer_data_string = "7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)"

    # 4. Parse the NMR data from the answer
    parsed_data = parse_nmr_data(answer_data_string)

    if parsed_data is None:
        return f"Incorrect. The NMR data for answer {answer_choice} ('{answer_data_string}') is malformed or could not be parsed."

    # 5. Check constraints
    # Check total proton integration
    total_integration = sum(s['integration'] for s in parsed_data)
    if total_integration != expected_total_protons:
        return f"Incorrect. The total number of protons is {total_integration}, but it should be {expected_total_protons} for a haloacetophenone (C8H7XO)."

    # Check for the methyl ketone signal (3H singlet)
    methyl_signals = [s for s in parsed_data if s['integration'] == 3 and s['multiplicity'] == 's']
    if not methyl_signals:
        return f"Incorrect. The spectrum is missing the characteristic 3H singlet for the acetyl group (-COCH3)."
    
    # Check the chemical shift of the methyl signal
    methyl_shift = methyl_signals[0]['shift']
    if not (2.0 <= methyl_shift <= 2.7):
        return f"Incorrect. The chemical shift of the 3H singlet ({methyl_shift} ppm) is outside the typical range (2.0-2.7 ppm) for a methyl ketone."

    # Check for aromatic signals
    aromatic_signals = [s for s in parsed_data if 6.5 <= s['shift'] <= 8.5]
    aromatic_proton_count = sum(s['integration'] for s in aromatic_signals)
    if aromatic_proton_count != expected_aromatic_protons:
        return f"Incorrect. The spectrum shows {aromatic_proton_count} aromatic protons, but a di-substituted benzene ring should have {expected_aromatic_protons}."

    # Check the specific pattern for para-substitution (two 2H doublets)
    aromatic_doublets_2h = [s for s in aromatic_signals if s['integration'] == 2 and s['multiplicity'] == 'd']
    if len(aromatic_doublets_2h) != 2:
        return f"Incorrect. The aromatic region does not show the characteristic pattern of two 2H doublets expected for a para-substituted ring. Found {len(aromatic_doublets_2h)} such signals."

    # If all checks pass
    return "Correct"

# Run the check
result = check_answer()
print(result)