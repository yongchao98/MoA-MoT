import re

def check_correctness():
    """
    This function checks if the provided answer 'C' (ethyl 4-aminobenzoate)
    is consistent with all the given spectroscopic data.
    """
    
    # --- Data from the Question ---
    GIVEN_FORMULA = "C9H11NO2"
    GIVEN_IR = {
        'amine_bands': 2,  # Two bands at 3420 and 3325 cm-1 imply a primary amine/amide
        'carbonyl_freq': 1720
    }
    GIVEN_NMR = {
        'signals': [
            {'ppm': 1.20, 'type': 't', 'H': 3},   # Triplet for -CH3
            {'ppm': 4.0,  'type': 'bs', 'H': 2},  # Broad singlet for -NH2
            {'ppm': 4.5,  'type': 'q', 'H': 2},   # Quartet for -CH2-
            {'ppm': 7.0,  'type': 'd', 'H': 2},   # Aromatic doublet
            {'ppm': 8.0,  'type': 'd', 'H': 2}    # Aromatic doublet
        ]
    }
    
    # --- Candidate Answer to Check ---
    # The final answer provided is <<<C>>>, which corresponds to ethyl 4-aminobenzoate.
    ANSWER_KEY = 'C'

    # --- Database of Candidate Properties ---
    # Note: Chemical shift and frequency ranges are typical approximations.
    compounds = {
        'A': {
            'name': '3-ethoxybenzamide',
            'formula': 'C9H11NO2',
            'amine_type': 'primary_amide', # -CONH2, gives 2 N-H bands
            'carbonyl_type': 'amide',
            'carbonyl_freq_range': (1650, 1690),
            'aromatic_pattern': 'meta', # 1,3-disubstituted -> complex pattern
            'ethyl_quartet_env': 'ether', # -O-CH2-Ar
            'ethyl_quartet_ppm_range': (3.9, 4.2)
        },
        'B': {
            'name': 'N-(4-ethoxyphenyl)formamide',
            'formula': 'C9H11NO2',
            'amine_type': 'secondary_amide', # -NH-, gives 1 N-H band
            'carbonyl_type': 'amide',
            'carbonyl_freq_range': (1660, 1700),
            'aromatic_pattern': 'para', # 1,4-disubstituted -> two doublets
            'ethyl_quartet_env': 'ether', # -O-CH2-Ar
            'ethyl_quartet_ppm_range': (3.9, 4.2)
        },
        'C': {
            'name': 'ethyl 4-aminobenzoate',
            'formula': 'C9H11NO2',
            'amine_type': 'primary_amine', # -NH2, gives 2 N-H bands
            'carbonyl_type': 'conjugated_ester',
            'carbonyl_freq_range': (1715, 1730),
            'aromatic_pattern': 'para', # 1,4-disubstituted -> two doublets
            'ethyl_quartet_env': 'ester_oxygen', # -COO-CH2-
            'ethyl_quartet_ppm_range': (4.2, 4.6)
        },
        'D': {
            'name': '4-aminophenyl propionate',
            'formula': 'C9H11NO2',
            'amine_type': 'primary_amine', # -NH2, gives 2 N-H bands
            'carbonyl_type': 'phenyl_ester',
            'carbonyl_freq_range': (1750, 1770), # Phenyl esters are higher
            'aromatic_pattern': 'para', # 1,4-disubstituted -> two doublets
            'ethyl_quartet_env': 'ester_carbonyl', # -CO-CH2-CH3
            'ethyl_quartet_ppm_range': (2.2, 2.6)
        }
    }

    chosen_compound = compounds[ANSWER_KEY]
    
    # 1. Check Molecular Formula
    if chosen_compound['formula'] != GIVEN_FORMULA:
        return f"Incorrect. The molecular formula of {chosen_compound['name']} ({chosen_compound['formula']}) does not match the given formula ({GIVEN_FORMULA})."

    # 2. Check IR - Amine/Amide Type
    if chosen_compound['amine_type'] in ['primary_amine', 'primary_amide'] and GIVEN_IR['amine_bands'] != 2:
        return f"Incorrect. The IR shows 2 N-H bands, but {chosen_compound['name']} is not a primary amine/amide."
    if chosen_compound['amine_type'] == 'secondary_amide' and GIVEN_IR['amine_bands'] != 1:
        return f"Incorrect. The IR shows 2 N-H bands (primary amine), but {chosen_compound['name']} is a secondary amide which would show only 1 band."

    # 3. Check IR - Carbonyl Frequency
    co_min, co_max = chosen_compound['carbonyl_freq_range']
    if not (co_min <= GIVEN_IR['carbonyl_freq'] <= co_max):
        return f"Incorrect. The IR carbonyl band is at {GIVEN_IR['carbonyl_freq']} cm-1, which is outside the expected range of {co_min}-{co_max} cm-1 for a {chosen_compound['carbonyl_type']} in {chosen_compound['name']}."

    # 4. Check NMR - Aromatic Pattern
    aromatic_signals = [s for s in GIVEN_NMR['signals'] if 6.5 <= s['ppm'] <= 8.5]
    is_para_pattern = len(aromatic_signals) == 2 and all(s['type'] == 'd' for s in aromatic_signals)
    if chosen_compound['aromatic_pattern'] == 'meta' and is_para_pattern:
        return f"Incorrect. The NMR shows a para-substitution pattern (two doublets), but {chosen_compound['name']} is meta-substituted and would show a more complex pattern."

    # 5. Check NMR - Ethyl Quartet Chemical Shift (Key Differentiator)
    quartet_signal_ppm = next(s['ppm'] for s in GIVEN_NMR['signals'] if s['type'] == 'q')
    q_min, q_max = chosen_compound['ethyl_quartet_ppm_range']
    if not (q_min <= quartet_signal_ppm <= q_max):
        return f"Incorrect. The NMR shows a quartet at {quartet_signal_ppm} ppm. This is characteristic of a -CH2- group attached to an ester oxygen. The expected shift for the -CH2- in {chosen_compound['name']} is {q_min}-{q_max} ppm, which does not match."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
print(check_correctness())