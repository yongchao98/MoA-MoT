import re

def check_nmr_answer():
    """
    This function programmatically checks the correctness of the LLM's answer
    by deducing the molecule's structure and its expected 1H NMR spectrum from
    the question's constraints.
    """
    # 1. Define constraints from the question and deduce the structure.
    # - 8 total carbons, di-substituted benzene ring -> C6H4 core + 2C in substituents.
    # - Carbonyl (C=O) and aromatic-halogen (-X) -> Substituents are -X and -COCH3.
    # - Deduced structure: A haloacetophenone (X-C6H4-COCH3).

    # 2. Define expected 1H NMR features for a haloacetophenone.
    expected_features = {
        "total_protons": 7,  # 4 on the ring + 3 on the methyl group
        "acetyl_methyl": {
            "integration": 3,
            "multiplicity": "s",
            "shift_range": (2.0, 2.7)
        },
        "aromatic_protons": {
            "integration": 4,
            "shift_range": (6.5, 8.5)
        },
        "disqualifying_signals": {
            "aldehyde_shift": 9.0 # A signal above this indicates an aldehyde, not a ketone.
        }
    }

    # 3. The LLM's answer and the corresponding data.
    llm_answer_key = "A"
    answer_data_string = "1H NMR: chemical reference (ppm): 7.8 (2H, d), 7.6 (2H, d), 2.3 (3H, s)"

    # 4. Parse the NMR data string into a structured format.
    def parse_nmr_string(s):
        signals = []
        # Regex to find patterns like "7.8 (2H, d)"
        pattern = re.compile(r"(\d+\.?\d*)\s*\((\d+)H,\s*([a-z]+)\)")
        matches = pattern.findall(s)
        for match in matches:
            signals.append({
                "shift": float(match[0]),
                "integration": int(match[1]),
                "multiplicity": match[2].strip()
            })
        return signals

    parsed_signals = parse_nmr_string(answer_data_string)

    # 5. Perform the verification checks.
    
    # Check for the acetyl methyl group signal
    acetyl_signal_found = False
    for signal in parsed_signals:
        if (signal["integration"] == expected_features["acetyl_methyl"]["integration"] and
            signal["multiplicity"] == expected_features["acetyl_methyl"]["multiplicity"] and
            expected_features["acetyl_methyl"]["shift_range"][0] <= signal["shift"] <= expected_features["acetyl_methyl"]["shift_range"][1]):
            acetyl_signal_found = True
            break
    if not acetyl_signal_found:
        return "Incorrect. The NMR data does not show a signal for an acetyl methyl group (a 3H singlet in the 2.0-2.7 ppm range)."

    # Check for the correct total number of aromatic protons
    aromatic_protons_total = 0
    for signal in parsed_signals:
        if expected_features["aromatic_protons"]["shift_range"][0] <= signal["shift"] <= expected_features["aromatic_protons"]["shift_range"][1]:
            aromatic_protons_total += signal["integration"]
    if aromatic_protons_total != expected_features["aromatic_protons"]["integration"]:
        return f"Incorrect. The total integration in the aromatic region is {aromatic_protons_total}H, but it should be {expected_features['aromatic_protons']['integration']}H for a di-substituted benzene ring."

    # Check for disqualifying signals (like an aldehyde)
    for signal in parsed_signals:
        if signal["shift"] >= expected_features["disqualifying_signals"]["aldehyde_shift"]:
            return f"Incorrect. A signal at {signal['shift']} ppm indicates an aldehyde, but the structure is a ketone."

    # Check total proton count
    total_protons_found = sum(s['integration'] for s in parsed_signals)
    if total_protons_found != expected_features["total_protons"]:
        return f"Incorrect. The total number of protons is {total_protons_found}, but the deduced structure (haloacetophenone) has {expected_features['total_protons']} protons."

    # Final check: The specific pattern of two doublets for 2H each is a hallmark of para-substitution,
    # which is a valid and common isomer for this structure. The data perfectly matches this.
    aromatic_signals = [s for s in parsed_signals if expected_features["aromatic_protons"]["shift_range"][0] <= s['shift']]
    if len(aromatic_signals) == 2 and all(s['multiplicity'] == 'd' and s['integration'] == 2 for s in aromatic_signals):
        # This confirms the highly symmetric para-substitution pattern.
        pass
    else:
        return "Incorrect. The aromatic signal pattern does not match the classic para-substitution pattern of two 2H doublets presented in the answer."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_nmr_answer()
print(result)