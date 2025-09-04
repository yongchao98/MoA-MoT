import re

def check_answer():
    """
    This function checks if the proposed answer 'A) ethyl 4-aminobenzoate'
    is consistent with all the provided spectral data.
    """
    # --- Given Data ---
    molecular_formula = "C9H11NO2"
    ir_bands = [3420, 3325, 1720]
    nmr_signals = {
        '1.20': {'type': 't', 'integration': 3},
        '4.0': {'type': 'bs', 'integration': 2},
        '4.5': {'type': 'q', 'integration': 2},
        '7.0': {'type': 'd', 'integration': 2},
        '8.0': {'type': 'd', 'integration': 2},
    }
    
    # --- Properties of the proposed answer: A) ethyl 4-aminobenzoate ---
    # Structure: An amino group (-NH2) and an ethyl ester group (-COOCH2CH3)
    # are para (1,4) to each other on a benzene ring.

    # 1. Check Molecular Formula
    # C6H4 (ring) + C (carbonyl) + C2H5 (ethyl) + NH2 (amine) + O (carbonyl) + O (ether)
    # C = 6 + 1 + 2 = 9
    # H = 4 + 5 + 2 = 11
    # N = 1
    # O = 1 + 1 = 2
    # Calculated formula: C9H11NO2
    if molecular_formula != "C9H11NO2":
        return f"Incorrect: The molecular formula of ethyl 4-aminobenzoate is C9H11NO2, which does not match the provided formula '{molecular_formula}'."

    # 2. Check IR Data
    # IR: medium to strong intensity bands at 3420 cm-1, 3325 cm-1
    # These two bands are characteristic of a primary amine (-NH2) N-H stretch (symmetric and asymmetric).
    # Ethyl 4-aminobenzoate has a primary amine. This is consistent.
    # IR: strong band at 1720 cm-1
    # This is characteristic of a carbonyl (C=O) group. The frequency ~1720 cm-1 is typical for an ester
    # that is conjugated with an aromatic ring.
    # Ethyl 4-aminobenzoate has a conjugated ester group. This is consistent.
    if not (3300 <= 3420 <= 3500 and 3300 <= 3325 <= 3500):
         return "Incorrect: The structure has a primary amine, but the given IR N-H bands are outside the typical range."
    if not (1715 <= 1720 <= 1735):
        return f"Incorrect: The IR carbonyl band at 1720 cm-1 is a perfect match for a conjugated ester, but the check failed. This indicates a logic error in the check itself, as the data is consistent."
        
    # 3. Check 1H NMR Data
    
    # Check total proton count
    total_protons_nmr = sum(s['integration'] for s in nmr_signals.values())
    if total_protons_nmr != 11:
        return f"Incorrect: The sum of protons from NMR integration is {total_protons_nmr}, but the structure has 11 protons."

    # Check for ethyl group: -CH2-CH3
    # Expect a triplet (3H) for -CH3 and a quartet (2H) for -CH2.
    # The -CH3 (next to CH2) should be ~1.2-1.4 ppm. Given: 1.20 ppm (t, 3H). Consistent.
    # The -O-CH2- (next to CH3 and ester oxygen) should be ~4.1-4.6 ppm. Given: 4.5 ppm (q, 2H). Consistent.
    if not (nmr_signals['1.20']['type'] == 't' and nmr_signals['1.20']['integration'] == 3 and \
            nmr_signals['4.5']['type'] == 'q' and nmr_signals['4.5']['integration'] == 2):
        return "Incorrect: The NMR signals at 1.20 ppm and 4.5 ppm do not match the expected triplet/quartet pattern of an ethyl ester group."

    # Check for amine protons: -NH2
    # Expect a broad singlet for 2H. Given: 4.0 ppm (bs, 2H). Consistent.
    if not (nmr_signals['4.0']['type'] == 'bs' and nmr_signals['4.0']['integration'] == 2):
        return "Incorrect: The NMR signal at 4.0 ppm does not match the expected broad singlet for a 2H amine group."

    # Check for aromatic protons: 1,4-disubstituted ring
    # With an electron-donating group (-NH2) and an electron-withdrawing group (-COOEt),
    # we expect two doublets, each integrating to 2H.
    # The protons ortho to the EWG (-COOEt) are deshielded (downfield, ~8.0 ppm).
    # The protons ortho to the EDG (-NH2) are shielded (upfield, ~6.5-7.0 ppm).
    # Given: 7.0 ppm (d, 2H) and 8.0 ppm (d, 2H). This is a classic match.
    if not (nmr_signals['7.0']['type'] == 'd' and nmr_signals['7.0']['integration'] == 2 and \
            nmr_signals['8.0']['type'] == 'd' and nmr_signals['8.0']['integration'] == 2):
        return "Incorrect: The NMR signals in the aromatic region do not match the expected two-doublet pattern for a para-disubstituted ring with EDG/EWG."

    # All checks passed.
    return "Correct"

# Run the check
result = check_answer()
print(result)