import collections

def check_correctness_of_answer():
    """
    This function checks if the proposed answer, ethyl 4-aminobenzoate,
    is consistent with all the provided spectral data (Molecular Formula, IR, 1H NMR).
    """
    # --- Data provided in the question ---
    given_data = {
        "formula": {'C': 9, 'H': 11, 'N': 1, 'O': 2},
        "ir_peaks": {
            "N-H_stretch": [3420, 3325],
            "C=O_stretch": 1720
        },
        "nmr_peaks": [
            {"shift": 1.20, "multiplicity": "t", "integration": 3},
            {"shift": 4.0, "multiplicity": "bs", "integration": 2},
            {"shift": 4.5, "multiplicity": "q", "integration": 2},
            {"shift": 7.0, "multiplicity": "d", "integration": 2},
            {"shift": 8.0, "multiplicity": "d", "integration": 2}
        ]
    }

    # --- Properties of the proposed answer: B) ethyl 4-aminobenzoate ---
    # Structure: A benzene ring with a primary amine (-NH2) and an ethyl ester (-COOCH2CH3)
    # at positions 1 and 4 (para).
    proposed_answer = {
        "name": "ethyl 4-aminobenzoate",
        "formula": {'C': 9, 'H': 11, 'N': 1, 'O': 2},
        "has_primary_amine": True,
        "carbonyl_type": "conjugated_ester",
        "nmr_features": {
            "has_ethyl_ester": True, # Expects -O-CH2-CH3
            "aromatic_pattern": "para_disubstituted"
        }
    }

    # 1. Check Molecular Formula
    if collections.Counter(proposed_answer["formula"]) != collections.Counter(given_data["formula"]):
        return f"Incorrect: The molecular formula of {proposed_answer['name']} ({proposed_answer['formula']}) does not match the given formula ({given_data['formula']})."

    # 2. Check IR Spectroscopy Data
    # Check for primary amine (-NH2)
    if proposed_answer["has_primary_amine"]:
        nh_stretches = given_data["ir_peaks"]["N-H_stretch"]
        if len(nh_stretches) != 2:
            return "Incorrect: The proposed structure is a primary amine, which requires two N-H stretch bands in the IR, but the data does not show two bands."
        for peak in nh_stretches:
            if not (3300 <= peak <= 3500):
                return f"Incorrect: The IR N-H stretch at {peak} cm-1 is outside the typical range (3300-3500 cm-1) for a primary amine."
    else:
        if "N-H_stretch" in given_data["ir_peaks"] and len(given_data["ir_peaks"]["N-H_stretch"]) == 2:
            return "Incorrect: The IR data shows two N-H bands, indicating a primary amine, but the proposed structure is not a primary amine."

    # Check for carbonyl (C=O) group
    co_stretch = given_data["ir_peaks"]["C=O_stretch"]
    if proposed_answer["carbonyl_type"] == "conjugated_ester":
        if not (1715 <= co_stretch <= 1730):
            return f"Incorrect: The IR C=O stretch at {co_stretch} cm-1 is not in the typical range for a conjugated ester (1715-1730 cm-1)."
    else:
        return f"Incorrect: The proposed structure's carbonyl type ({proposed_answer['carbonyl_type']}) is not consistent with the IR peak at {co_stretch} cm-1."

    # 3. Check 1H NMR Spectroscopy Data
    nmr_peaks = given_data["nmr_peaks"]
    
    # Check total integration
    total_integration = sum(p['integration'] for p in nmr_peaks)
    if total_integration != given_data['formula']['H']:
        return f"Incorrect: The total integration of the 1H NMR spectrum ({total_integration}H) does not match the number of hydrogens in the molecular formula ({given_data['formula']['H']}H)."

    # Check for ethyl ester group (-O-CH2-CH3)
    triplet_3h = any(p['shift'] > 1.1 and p['shift'] < 1.4 and p['multiplicity'] == 't' and p['integration'] == 3 for p in nmr_peaks)
    quartet_2h_ester = any(p['shift'] > 4.2 and p['shift'] < 4.6 and p['multiplicity'] == 'q' and p['integration'] == 2 for p in nmr_peaks)
    if not (triplet_3h and quartet_2h_ester):
        return "Incorrect: The 1H NMR data lacks the characteristic signals for an ethyl ester group: a triplet (3H) around 1.2-1.3 ppm and a quartet (2H) around 4.3-4.5 ppm."

    # Check for primary amine (-NH2)
    amine_signal = any(p['shift'] > 3.0 and p['shift'] < 5.0 and p['multiplicity'] in ['bs', 's'] and p['integration'] == 2 for p in nmr_peaks)
    if not amine_signal:
        return "Incorrect: The 1H NMR data lacks a signal consistent with a primary amine (typically a broad singlet for 2H)."

    # Check for para-disubstituted aromatic ring
    aromatic_signals = [p for p in nmr_peaks if p['shift'] > 6.5 and p['shift'] < 8.5]
    if len(aromatic_signals) == 2:
        is_para_pattern = all(p['multiplicity'] == 'd' and p['integration'] == 2 for p in aromatic_signals)
        if not is_para_pattern:
            return "Incorrect: The aromatic region of the 1H NMR spectrum is not consistent with a 1,4-disubstituted (para) benzene ring, which should show two doublets of 2H each."
    else:
        return "Incorrect: The number of signals in the aromatic region of the 1H NMR spectrum is not consistent with a 1,4-disubstituted (para) benzene ring."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness_of_answer()
print(result)