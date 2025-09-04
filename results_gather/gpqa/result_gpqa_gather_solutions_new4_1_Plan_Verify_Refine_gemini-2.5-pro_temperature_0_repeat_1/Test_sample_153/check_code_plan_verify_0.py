import sys
from io import StringIO

def check_answer():
    """
    This function checks the correctness of the proposed answer by analyzing the given spectral data.
    """
    # --- Define Spectral Data and Options ---
    question_data = {
        "ms": {"m_z_156": 100, "m_z_158": 32},
        "ir": {"broad_oh": (3500, 2700), "carbonyl": 1720},
        "nmr": {
            "signals": [
                {"ppm": 11.0, "splitting": "s", "integration": 1, "type": "acidic"},
                {"ppm": 8.02, "splitting": "d", "integration": 2, "type": "aromatic"},
                {"ppm": 7.72, "splitting": "d", "integration": 2, "type": "aromatic"}
            ]
        }
    }

    # Define properties of the possible answers
    compounds = {
        'A': {'name': '3-Chloro-2-hydroxybenzaldehyde', 'formula': {'C': 7, 'H': 5, 'Cl': 1, 'O': 2}, 'functional_groups': ['phenol', 'aldehyde'], 'substitution': '1,2,3-trisubstituted', 'aromatic_protons': 3, 'acidic_protons': 1},
        'B': {'name': 'Phenyl chloroformate', 'formula': {'C': 7, 'H': 5, 'Cl': 1, 'O': 2}, 'functional_groups': ['chloroformate', 'ester'], 'substitution': 'monosubstituted', 'aromatic_protons': 5, 'acidic_protons': 0},
        'C': {'name': '2-chlorobenzoic acid', 'formula': {'C': 7, 'H': 5, 'Cl': 1, 'O': 2}, 'functional_groups': ['carboxylic_acid'], 'substitution': 'ortho', 'aromatic_protons': 4, 'acidic_protons': 1},
        'D': {'name': '4-chlorobenzoic acid', 'formula': {'C': 7, 'H': 5, 'Cl': 1, 'O': 2}, 'functional_groups': ['carboxylic_acid'], 'substitution': 'para', 'aromatic_protons': 4, 'acidic_protons': 1}
    }

    # The final answer provided by the LLM
    final_answer_key = 'D'
    proposed_compound = compounds[final_answer_key]
    
    # --- Verification Checks ---
    
    # 1. Mass Spectrometry Check
    # Check for one chlorine atom
    ms_ratio = question_data["ms"]["m_z_158"] / question_data["ms"]["m_z_156"]
    # Expected ratio for one Cl is ~0.33 (24.2% / 75.8%)
    if not (0.30 <= ms_ratio <= 0.35):
        return f"Incorrect. The MS M+2/M+ ratio of {ms_ratio:.2f} does not match the expected ~0.33 for a single chlorine atom."
    if proposed_compound['formula'].get('Cl', 0) != 1:
        return f"Incorrect. The MS data indicates one chlorine atom, but the proposed structure '{proposed_compound['name']}' does not contain exactly one."
        
    # Check molecular weight for the M+ peak (using 35Cl)
    atomic_weights = {'C': 12, 'H': 1, 'O': 16, 'Cl': 35}
    formula = proposed_compound['formula']
    mw = sum(atomic_weights[atom] * count for atom, count in formula.items())
    if mw != 156:
        return f"Incorrect. The MS M+ peak is at m/z=156, but the calculated molecular weight for '{proposed_compound['name']}' (using 35Cl) is {mw}."

    # 2. IR Spectroscopy Check
    # Check for carboxylic acid features
    if 'carboxylic_acid' not in proposed_compound['functional_groups']:
        return f"Incorrect. The IR spectrum (broad peak at 3500-2700 cm-1 and sharp peak at 1720 cm-1) is characteristic of a carboxylic acid. The proposed structure, '{proposed_compound['name']}', is not a carboxylic acid."

    # 3. 1H NMR Spectroscopy Check
    # Check for carboxylic acid proton
    has_acid_proton_signal = any(s['type'] == 'acidic' and 10 <= s['ppm'] <= 13 for s in question_data['nmr']['signals'])
    if has_acid_proton_signal and 'carboxylic_acid' not in proposed_compound['functional_groups']:
        return f"Incorrect. The 1H NMR spectrum shows a signal at 11.0 ppm, characteristic of a carboxylic acid proton. The proposed structure, '{proposed_compound['name']}', does not have this proton."

    # Check aromatic proton count
    total_aromatic_protons_data = sum(s['integration'] for s in question_data['nmr']['signals'] if s['type'] == 'aromatic')
    if total_aromatic_protons_data != proposed_compound['aromatic_protons']:
        return f"Incorrect. The 1H NMR spectrum shows a total of {total_aromatic_protons_data} aromatic protons. The proposed structure, '{proposed_compound['name']}', has {proposed_compound['aromatic_protons']} aromatic protons."

    # Check aromatic substitution pattern
    aromatic_signals = [s for s in question_data['nmr']['signals'] if s['type'] == 'aromatic']
    # The pattern is two doublets, each for 2H. This is a classic para-substitution pattern.
    is_para_pattern = (
        len(aromatic_signals) == 2 and
        all(s['splitting'] == 'd' for s in aromatic_signals) and
        all(s['integration'] == 2 for s in aromatic_signals)
    )
    if is_para_pattern and proposed_compound['substitution'] != 'para':
        return f"Incorrect. The 1H NMR spectrum shows a pattern of two doublets (2H each), which is characteristic of a para-substituted benzene ring. The proposed structure, '{proposed_compound['name']}', is {proposed_compound['substitution']}-substituted and would produce a different splitting pattern."

    # If all checks pass
    return "Correct"

# Run the check and print the result
result = check_answer()
print(result)