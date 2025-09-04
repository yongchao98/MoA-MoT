import re

def check_spectroscopy_answer():
    """
    This function checks if the proposed structure, ethyl 4-aminobenzoate,
    is consistent with all the provided spectral data.
    """
    # --- Given Data from the Question ---
    given_formula_str = "C9H11NO2"
    given_ir = {
        "N-H_stretches": [3420, 3325],
        "C=O_stretch": 1720
    }
    given_nmr = [
        {'ppm': 1.20, 'splitting': 't', 'integration': 3, 'label': 'CH3'},
        {'ppm': 4.0, 'splitting': 'bs', 'integration': 2, 'label': 'NH2'},
        {'ppm': 4.5, 'splitting': 'q', 'integration': 2, 'label': 'CH2'},
        {'ppm': 7.0, 'splitting': 'd', 'integration': 2, 'label': 'Aromatic H'},
        {'ppm': 8.0, 'splitting': 'd', 'integration': 2, 'label': 'Aromatic H'},
    ]
    
    # --- Properties of the Proposed Answer (D) ethyl 4-aminobenzoate ---
    proposed_structure = {
        "name": "ethyl 4-aminobenzoate",
        "formula": {'C': 9, 'H': 11, 'N': 1, 'O': 2},
        "functional_groups": ["primary aromatic amine", "conjugated ester", "ethyl group", "para-disubstituted benzene"]
    }

    # 1. Check Molecular Formula and Degree of Unsaturation (DoU)
    try:
        match = re.match(r"C(\d+)H(\d+)N(\d+)O(\d+)", given_formula_str)
        c, h, n, o = map(int, match.groups())
        if proposed_structure["formula"] != {'C': c, 'H': h, 'N': n, 'O': o}:
            return f"Incorrect Formula: The formula of {proposed_structure['name']} ({proposed_structure['formula']}) does not match the given formula {given_formula_str}."
        
        dou = c + 1 - (h / 2) + (n / 2)
        # Expected DoU for ethyl 4-aminobenzoate: 4 (benzene ring) + 1 (C=O) = 5
        if dou != 5:
            return f"Incorrect DoU: The Degree of Unsaturation for the formula {given_formula_str} is {dou}, which is inconsistent with the analysis in the answer (should be 5)."
    except (AttributeError, ValueError):
        return "Invalid format for the given molecular formula string."

    # 2. Check IR Spectroscopy Data
    # Check for primary amine (-NH2): two bands in the 3300-3500 cm-1 region.
    if "primary aromatic amine" not in proposed_structure["functional_groups"]:
        return "IR Mismatch: The proposed structure lacks a primary amine, but the IR shows two N-H bands at 3420 and 3325 cm-1."
    if len(given_ir["N-H_stretches"]) != 2:
        return "IR Mismatch: A primary amine should show two N-H bands, but the data describes a different number."

    # Check for conjugated ester (C=O): strong band ~1710-1730 cm-1.
    if "conjugated ester" not in proposed_structure["functional_groups"]:
        return "IR Mismatch: The proposed structure lacks a conjugated ester, but the IR shows a strong C=O band at 1720 cm-1."
    if not (1710 <= given_ir["C=O_stretch"] <= 1735):
        return f"IR Mismatch: The C=O frequency of {given_ir['C=O_stretch']} cm-1 is outside the typical range for a conjugated ester."

    # 3. Check 1H NMR Spectroscopy Data
    # Check total integration vs formula
    total_integration = sum(s['integration'] for s in given_nmr)
    if total_integration != proposed_structure["formula"]['H']:
        return f"NMR Mismatch: Total proton integration ({total_integration}H) does not match the number of hydrogens in the formula ({proposed_structure['formula']['H']}H)."

    # Check for ethyl group (-CH2CH3) attached to an oxygen
    triplet_3h = next((s for s in given_nmr if s['splitting'] == 't' and s['integration'] == 3 and 1.1 < s['ppm'] < 1.5), None)
    quartet_2h = next((s for s in given_nmr if s['splitting'] == 'q' and s['integration'] == 2 and 4.0 < s['ppm'] < 4.7), None)
    if not (triplet_3h and quartet_2h):
        return "NMR Mismatch: The spectrum does not contain the characteristic signals for an ethyl group attached to an oxygen (a triplet ~1.2 ppm and a quartet ~4.5 ppm)."

    # Check for primary amine (-NH2)
    amine_signal = next((s for s in given_nmr if s['splitting'] == 'bs' and s['integration'] == 2), None)
    if not amine_signal:
        return "NMR Mismatch: The spectrum lacks a broad singlet for 2H, which is expected for the primary amine group (-NH2)."

    # Check for para-disubstituted benzene ring
    aromatic_signals = [s for s in given_nmr if 6.5 < s['ppm'] < 8.5]
    if len(aromatic_signals) != 2:
        return f"NMR Mismatch: Expected 2 signals in the aromatic region for a para-disubstituted ring, but found {len(aromatic_signals)}."
    
    doublets_2h = [s for s in aromatic_signals if s['splitting'] == 'd' and s['integration'] == 2]
    if len(doublets_2h) != 2:
        return "NMR Mismatch: The aromatic signals are not two doublets of 2H each, which is the classic pattern for a 1,4-disubstituted benzene ring with an EDG and EWG."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_spectroscopy_answer()
print(result)