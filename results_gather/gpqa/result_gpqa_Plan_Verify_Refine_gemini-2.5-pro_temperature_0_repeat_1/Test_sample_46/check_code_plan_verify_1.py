import collections

def check_answer():
    """
    This function checks the correctness of the identified compound based on the provided spectral data.
    """
    # --- Given Data from the Question ---
    molecular_formula = "C9H11NO2"
    ir_data = {
        'N-H': [3420, 3325],  # Two bands indicate primary amine/amide
        'C=O': 1720
    }
    nmr_data = [
        {'ppm': 1.20, 'multiplicity': 't', 'integration': 3},
        {'ppm': 4.0, 'multiplicity': 'bs', 'integration': 2},
        {'ppm': 4.5, 'multiplicity': 'q', 'integration': 2},
        {'ppm': 7.0, 'multiplicity': 'd', 'integration': 2},
        {'ppm': 8.0, 'multiplicity': 'd', 'integration': 2},
    ]
    
    # --- Proposed Answer ---
    proposed_answer = 'B' # ethyl 4-aminobenzoate

    # --- Define Properties of Each Option ---
    # Chemical shift ranges are approximate.
    structures = {
        'A': {
            'name': '3-ethoxybenzamide',
            'formula': 'C9H11NO2',
            'ir_nh_bands': 2,  # Primary amide
            'ir_co_range': (1660, 1690), # Conjugated amide
            'nmr_aromatic_pattern': 'meta', # 1,3-disubstituted, complex pattern expected
            'nmr_ethyl_quartet_range': (3.9, 4.2) # -OCH2-
        },
        'B': {
            'name': 'ethyl 4-aminobenzoate',
            'formula': 'C9H11NO2',
            'ir_nh_bands': 2, # Primary amine
            'ir_co_range': (1715, 1730), # Conjugated ester
            'nmr_aromatic_pattern': 'para', # 1,4-disubstituted, two doublets expected
            'nmr_ethyl_quartet_range': (4.2, 4.6) # -COOCH2-
        },
        'C': {
            'name': 'N-(4-ethoxyphenyl)formamide',
            'formula': 'C9H11NO2',
            'ir_nh_bands': 1, # Secondary amide
            'ir_co_range': (1670, 1700),
            'nmr_aromatic_pattern': 'para',
            'nmr_has_formyl_H': True # Expects a signal > 8.0 ppm for -CHO
        },
        'D': {
            'name': '4-aminophenyl propionate',
            'formula': 'C9H11NO2',
            'ir_nh_bands': 2, # Primary amine
            'ir_co_range': (1735, 1755), # Phenyl ester
            'nmr_aromatic_pattern': 'para',
            'nmr_alkyl_quartet_range': (2.4, 2.8) # -COCH2-
        }
    }

    # --- Verification Logic ---
    
    # Get the properties of the proposed answer
    props = structures[proposed_answer]

    # 1. Check Molecular Formula
    # All options have the correct formula C9H11NO2, so this check is trivial.
    # Let's calculate it for the proposed answer to be sure.
    # ethyl 4-aminobenzoate: C6H4(NH2)(COOC2H5) -> C(6+1+2) H(4+2+5) N(1) O(2) -> C9H11NO2
    if props['formula'] != molecular_formula:
        return f"Incorrect. The molecular formula of {props['name']} is not {molecular_formula}."

    # 2. Check IR Data
    # Check N-H stretches
    if len(ir_data['N-H']) != props['ir_nh_bands']:
        return f"Incorrect. The IR data shows {len(ir_data['N-H'])} N-H bands, but {props['name']} (a { 'primary' if props['ir_nh_bands']==2 else 'secondary' } amine/amide) should have {props['ir_nh_bands']} band(s)."
    
    # Check C=O stretch
    co_min, co_max = props['ir_co_range']
    if not (co_min <= ir_data['C=O'] <= co_max):
        return f"Incorrect. The IR C=O band at {ir_data['C=O']} cm-1 is outside the expected range of {co_min}-{co_max} cm-1 for {props['name']}."

    # 3. Check 1H NMR Data
    # Check for ethyl group pattern: 3H triplet and 2H quartet
    signals = sorted(nmr_data, key=lambda x: x['ppm'])
    triplet_3h = any(s['multiplicity'] == 't' and s['integration'] == 3 for s in signals)
    quartet_2h = any(s['multiplicity'] == 'q' and s['integration'] == 2 for s in signals)
    if not (triplet_3h and quartet_2h):
        return "Incorrect. The NMR data does not show the characteristic triplet-quartet pattern of an ethyl group."
    
    # Check the chemical shift of the quartet
    quartet_signal = next(s for s in signals if s['multiplicity'] == 'q' and s['integration'] == 2)
    q_min, q_max = props['nmr_ethyl_quartet_range']
    if not (q_min <= quartet_signal['ppm'] <= q_max):
        # Let's check why other options would be wrong here
        if proposed_answer != 'D' and (2.4 <= quartet_signal['ppm'] <= 2.8):
             return f"Incorrect. The quartet at {quartet_signal['ppm']} ppm suggests a -COCH2- group (like in option D), not the -COOCH2- group expected for {props['name']}."
        return f"Incorrect. The quartet at {quartet_signal['ppm']} ppm is outside the expected range of {q_min}-{q_max} ppm for the -CH2- group in {props['name']}."

    # Check for aromatic pattern
    aromatic_doublets = [s for s in signals if s['multiplicity'] == 'd' and s['integration'] == 2 and s['ppm'] > 6.5]
    if props['nmr_aromatic_pattern'] == 'para' and len(aromatic_doublets) != 2:
        return f"Incorrect. The NMR data does not show two doublets in the aromatic region, which is expected for a 1,4-disubstituted (para) ring like in {props['name']}."
    if props['nmr_aromatic_pattern'] == 'meta' and len(aromatic_doublets) == 2:
        return f"Incorrect. The NMR data shows two clean doublets, which indicates a para-substituted ring, not a meta-substituted ring as in {props['name']}."

    # Check for amine protons
    amine_signal = any(s['multiplicity'] == 'bs' and s['integration'] == 2 for s in signals)
    if not amine_signal:
        return f"Incorrect. The NMR data lacks a broad singlet for 2H, which is expected for the primary amine (-NH2) group in {props['name']}."

    # All checks for the proposed answer passed.
    return "Correct"

# Run the check
result = check_answer()
print(result)