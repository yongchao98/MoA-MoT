import collections

def find_room_temp_lcs():
    """
    Analyzes experimental data for n-alkyl cyanobenzene (nCN) liquid crystals
    to find candidates that are in a liquid crystal phase at room temperature.
    """
    target_temp = 25.0  # Room temperature in Celsius

    # Experimental data for n-alkyl cyanobenzenes (CnH2n+1-Ph-CN)
    # Data includes melting point (Cr -> Sm/N) and clearing points (Sm -> N, N -> I)
    # 'phases' is a list of (temp, phase_from, phase_to) tuples
    lc_data = [
        {
            'n': 5,
            'name': '4-pentylbenzonitrile',
            'formula': 'C5H11-Ph-CN',
            'transitions': [(-7, 'Cr', 'N'), (-2, 'N', 'I')]
        },
        {
            'n': 6,
            'name': '4-hexylbenzonitrile',
            'formula': 'C6H13-Ph-CN',
            'transitions': [(13.9, 'Cr', 'N'), (29.5, 'N', 'I')]
        },
        {
            'n': 7,
            'name': '4-heptylbenzonitrile',
            'formula': 'C7H15-Ph-CN',
            'transitions': [(30.2, 'Cr', 'N'), (42.8, 'N', 'I')]
        },
        {
            'n': 8,
            'name': '4-octylbenzonitrile',
            'formula': 'C8H17-Ph-CN',
            'transitions': [(21.5, 'Cr', 'SmA'), (33.5, 'SmA', 'N'), (40.5, 'N', 'I')]
        },
        {
            'n': 9,
            'name': '4-nonylbenzonitrile',
            'formula': 'C9H19-Ph-CN',
            'transitions': [(42.0, 'Cr', 'SmA'), (47.5, 'SmA', 'N'), (49.5, 'N', 'I')]
        }
    ]

    print(f"Searching for single-ring liquid crystals with a phase at {target_temp}°C...\n")
    
    found_candidates = []

    for molecule in lc_data:
        # Sort transitions by temperature to ensure correct phase ranges
        sorted_transitions = sorted(molecule['transitions'])
        
        for i in range(len(sorted_transitions) - 1):
            temp_low, phase_low_from, phase_low_to = sorted_transitions[i]
            temp_high, _, _ = sorted_transitions[i+1]
            
            # Check if the target temperature is within a liquid crystal phase range
            if temp_low <= target_temp < temp_high:
                found_candidates.append({
                    'name': molecule['name'],
                    'formula': molecule['formula'],
                    'n': molecule['n'],
                    'phase_at_target': phase_low_to,
                    'range': f"{temp_low}°C to {temp_high}°C",
                    'full_sequence': ' -> '.join([f"{t}°C ({p_to})" for t, p_from, p_to in sorted_transitions])
                })

    if found_candidates:
        print("Found suitable candidates:\n")
        for candidate in found_candidates:
            print(f"  Molecule: {candidate['name']} (n={candidate['n']})")
            print(f"  Formula: {candidate['formula']}")
            print(f"  At {target_temp}°C, the phase is: {candidate['phase_at_target']}")
            print(f"  This phase is stable in the range: {candidate['range']}")
            print(f"  Full phase sequence: Cr -> {candidate['full_sequence']}")
            print("-" * 30)
    else:
        print("No suitable candidates found in the dataset for the target temperature.")


if __name__ == '__main__':
    find_room_temp_lcs()
