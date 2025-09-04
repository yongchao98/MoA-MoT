import re

def check_spectroscopy_answer():
    """
    Checks if the proposed structure, ethyl 4-aminobenzoate, matches the given spectral data.
    """
    # --- 1. Define Given Data and Proposed Structure ---
    
    # Data from the question
    molecular_formula = "C9H11NO2"
    ir_peaks = {'N-H': [3420, 3325], 'C=O': 1720}
    nmr_signals = {
        1.20: {'mult': 't', 'int': 3, 'desc': 'Alkyl CH3'},
        4.0: {'mult': 'bs', 'int': 2, 'desc': 'Amine NH2'},
        4.5: {'mult': 'q', 'int': 2, 'desc': 'Alkyl CH2'},
        7.0: {'mult': 'd', 'int': 2, 'desc': 'Aromatic CH'},
        8.0: {'mult': 'd', 'int': 2, 'desc': 'Aromatic CH'}
    }
    
    # Properties of the proposed answer: B) ethyl 4-aminobenzoate
    structure = {
        'name': 'ethyl 4-aminobenzoate',
        'formula': 'C9H11NO2',
        'features': {
            'amine_type': 'primary', # -NH2 group
            'carbonyl_type': 'aryl ester', # -COO-Ar
            'substitution': 'para', # 1,4-disubstituted ring
            'fragments': ['ethyl_group'] # -CH2CH3
        }
    }

    # --- 2. Perform Checks ---

    # Check 1: Molecular Formula
    if structure['formula'] != molecular_formula:
        return f"Incorrect: The molecular formula of {structure['name']} is {structure['formula']}, which does not match the given formula {molecular_formula}."

    # Check 2: Degree of Unsaturation (DoU)
    # DoU = C + 1 - H/2 + N/2
    atoms = re.findall(r'([A-Z])(\d+)', molecular_formula)
    C, H, N = int(atoms[0][1]), int(atoms[1][1]), int(atoms[2][1])
    dou = C + 1 - (H / 2) + (N / 2)
    # Expected for ethyl 4-aminobenzoate: 4 (benzene ring) + 1 (C=O) = 5
    if dou != 5:
        return f"Incorrect: The calculated Degree of Unsaturation is {dou}, which is inconsistent with a benzene ring and a carbonyl group."

    # Check 3: IR Spectroscopy Constraints
    # Constraint 3a: N-H stretch. Two bands at 3420/3325 cm-1 strongly indicate a primary amine (-NH2).
    if structure['features']['amine_type'] != 'primary':
        return f"Incorrect: The IR spectrum shows two N-H bands (3420, 3325 cm-1), characteristic of a primary amine. The proposed structure is not a primary amine."
    
    # Constraint 3b: C=O stretch. A strong band at 1720 cm-1.
    # This is too high for an amide (~1650-1690 cm-1) and fits perfectly for a conjugated/aryl ester (~1715-1730 cm-1).
    if structure['features']['carbonyl_type'] != 'aryl ester':
        return f"Incorrect: The IR C=O band at 1720 cm-1 is characteristic of an aryl ester, not an amide or other carbonyl type."

    # Check 4: 1H NMR Spectroscopy Constraints
    # Constraint 4a: Aromatic region. Two doublets (7.0, 8.0 ppm), each for 2H.
    # This is a classic pattern for a 1,4- (para) disubstituted benzene ring.
    aromatic_doublets = [s for s in nmr_signals.values() if s['mult'] == 'd' and s['int'] == 2]
    if len(aromatic_doublets) != 2 or structure['features']['substitution'] != 'para':
        return f"Incorrect: The 1H NMR shows two doublets for 2H each in the aromatic region, indicating a para-substituted ring. The proposed structure does not match this pattern."

    # Constraint 4b: Ethyl group. A triplet (3H) and a quartet (2H).
    # Data: 1.20 ppm (t, 3H) and 4.5 ppm (q, 2H).
    triplet_found = any(s['mult'] == 't' and s['int'] == 3 for s in nmr_signals.values())
    quartet_found = any(s['mult'] == 'q' and s['int'] == 2 for s in nmr_signals.values())
    if not (triplet_found and quartet_found and 'ethyl_group' in structure['features']['fragments']):
        return f"Incorrect: The NMR data clearly shows an ethyl group (triplet and quartet). The proposed structure must contain this fragment."
    
    # Constraint 4c: Chemical shift of the ethyl quartet.
    # The quartet is at a very downfield 4.5 ppm, indicating it's attached to an oxygen (-O-CH2-), as in an ester.
    # A propionyl group (-CO-CH2-) would be around 2.5 ppm.
    quartet_shift = [k for k, v in nmr_signals.items() if v['mult'] == 'q'][0]
    if not (4.1 <= quartet_shift <= 4.6):
         return f"Incorrect: The quartet at {quartet_shift} ppm is too far downfield for a group attached to a carbonyl (like in a propionate) and is characteristic of a group attached to an ester oxygen (-O-CH2-)."

    # Constraint 4d: Amine protons. A broad singlet at 4.0 ppm for 2H.
    # This is consistent with the -NH2 protons of a primary aromatic amine.
    amine_signal_found = any(s['mult'] == 'bs' and s['int'] == 2 for s in nmr_signals.values())
    if not amine_signal_found:
        return f"Incorrect: The NMR data shows a broad singlet for 2H, consistent with amine protons, which is missing from the proposed structure's expected spectrum."

    # --- 3. Final Verdict ---
    # If all checks pass, the structure is consistent with all data.
    return "Correct"

# Execute the check
result = check_spectroscopy_answer()
print(result)