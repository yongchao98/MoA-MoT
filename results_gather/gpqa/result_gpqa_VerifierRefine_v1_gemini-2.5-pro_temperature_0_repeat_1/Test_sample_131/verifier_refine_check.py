import math
from collections import Counter

def check_correctness():
    """
    Checks the correctness of the answer by simulating the 1H NMR spectrum
    of the proposed mixture and comparing it to the given data.
    """

    # Step 1: Define the predicted 1H NMR data for each possible compound.
    # This data is derived from the chemical structure and symmetry of each molecule.
    compounds_db = {
        '1,2,4,5-tetramethylbenzene': {
            'name': '1,2,4,5-tetramethylbenzene (Durene)',
            'formula': 'C10H14',
            'aromatic_signals': [{'multiplicity': 'singlet', 'integration': 2}], # 2 equivalent aromatic H
            'alkyl_signals': [{'multiplicity': 'singlet', 'integration': 12}] # 4 equivalent methyl groups
        },
        '1,2,3,5-tetramethylbenzene': {
            'name': '1,2,3,5-tetramethylbenzene (Isodurene)',
            'formula': 'C10H14',
            'aromatic_signals': [{'multiplicity': 'singlet', 'integration': 2}], # 2 equivalent aromatic H
            'alkyl_signals': [
                {'multiplicity': 'singlet', 'integration': 6}, # 2 equivalent methyl groups
                {'multiplicity': 'singlet', 'integration': 3}, # 1 unique methyl group
                {'multiplicity': 'singlet', 'integration': 3}  # 1 unique methyl group
            ]
        },
        '1,2,3,4-tetramethylbenzene': {
            'name': '1,2,3,4-tetramethylbenzene (Prehnitene)',
            'formula': 'C10H14',
            'aromatic_signals': [{'multiplicity': 'singlet', 'integration': 2}], # 2 equivalent aromatic H
            'alkyl_signals': [
                {'multiplicity': 'singlet', 'integration': 6}, # 2 equivalent methyl groups
                {'multiplicity': 'singlet', 'integration': 6}  # 2 other equivalent methyl groups
            ]
        },
        '1,4-diethylbenzene': {
            'name': '1,4-diethylbenzene',
            'formula': 'C10H14',
            'aromatic_signals': [{'multiplicity': 'singlet', 'integration': 4}], # 4 equivalent aromatic H
            'alkyl_signals': [
                {'multiplicity': 'quartet', 'integration': 4}, # -CH2- groups
                {'multiplicity': 'triplet', 'integration': 6}  # -CH3 groups
            ]
        }
    }

    # The proposed answer is D, a mixture of 1,2,4,5-tetramethylbenzene and 1,2,3,4-tetramethylbenzene.
    compound1_name = '1,2,4,5-tetramethylbenzene'
    compound2_name = '1,2,3,4-tetramethylbenzene'

    c1 = compounds_db[compound1_name]
    c2 = compounds_db[compound2_name]

    # Step 2: Verify the mixture against the constraints from the question.

    # --- Aromatic Region Check (2 singlets, 1:1 ratio) ---
    aromatic_signals = c1['aromatic_signals'] + c2['aromatic_signals']

    if len(aromatic_signals) != 2:
        return f"Incorrect: The mixture of {c1['name']} and {c2['name']} produces {len(aromatic_signals)} aromatic signals, but the spectrum requires 2."
    
    if not all(s['multiplicity'] == 'singlet' for s in aromatic_signals):
        return f"Incorrect: Not all aromatic signals in the proposed mixture are singlets."

    aromatic_integrations = [s['integration'] for s in aromatic_signals]
    if aromatic_integrations.count(aromatic_integrations[0]) != len(aromatic_integrations):
        # This checks if all integrations are equal, which is required for a 1:1 ratio.
        return f"Incorrect: The aromatic signals have integrations {aromatic_integrations}, which is not a 1:1 ratio."

    # --- Alkyl Region Check (3 singlets, 2:1:1 ratio) ---
    alkyl_signals = c1['alkyl_signals'] + c2['alkyl_signals']

    if not all(s['multiplicity'] == 'singlet' for s in alkyl_signals):
        # This is a crucial check that eliminates 1,4-diethylbenzene.
        return f"Incorrect: Not all alkyl signals in the proposed mixture are singlets. The spectrum shows only singlets."

    if len(alkyl_signals) != 3:
        # This check assumes no accidental overlap of signals, which is standard for this type of problem.
        # For the correct answer, 1,2,4,5-TMB has 1 signal and 1,2,3,4-TMB has 2 signals, for a total of 3.
        return f"Incorrect: The mixture produces {len(alkyl_signals)} distinct alkyl signals, but the spectrum requires 3."

    # Check the 2:1:1 ratio. For a 1:1 mixture, the total alkyl protons are 12+12=24.
    # A 2:1:1 ratio corresponds to integrations of 12H, 6H, and 6H.
    alkyl_integrations = sorted([s['integration'] for s in alkyl_signals], reverse=True)
    expected_integrations = [12, 6, 6]

    if alkyl_integrations != expected_integrations:
        # We can also check the ratio directly
        common_divisor = math.gcd(math.gcd(alkyl_integrations[0], alkyl_integrations[1]), alkyl_integrations[2])
        ratio = sorted([i // common_divisor for i in alkyl_integrations], reverse=True)
        return f"Incorrect: The alkyl signals have an integration ratio of {':'.join(map(str, ratio))}, not 2:1:1."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)