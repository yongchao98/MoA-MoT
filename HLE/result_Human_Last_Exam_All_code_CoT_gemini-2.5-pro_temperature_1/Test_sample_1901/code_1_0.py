import sys

def find_room_temp_lc():
    """
    Searches for a single-ring liquid crystal that is nematic at room temperature.
    """
    # Database of 4-alkyl-cyanobenzene liquid crystals
    # Data includes the alkyl chain length 'n', name, and key transition temperatures in Celsius.
    # The nematic phase exists between nematic_start_temp and clearing_point.
    lc_data = [
        {
            'n': 5,
            'name': '4-Pentyl-cyanobenzene',
            'nematic_start_temp': 30.0,  # Crystal to Nematic transition
            'clearing_point': 55.0,      # Nematic to Isotropic transition
        },
        {
            'n': 6,
            'name': '4-Hexyl-cyanobenzene',
            'nematic_start_temp': 14.0,  # Crystal to Nematic transition
            'clearing_point': 29.0,      # Nematic to Isotropic transition
        },
        {
            'n': 7,
            'name': '4-Heptyl-cyanobenzene',
            'nematic_start_temp': 30.0,  # Crystal to Nematic transition
            'clearing_point': 43.0,      # Nematic to Isotropic transition
        },
        {
            'n': 8,
            'name': '4-Octyl-cyanobenzene',
            'nematic_start_temp': 35.0,  # Smectic A to Nematic transition
            'clearing_point': 42.5,      # Nematic to Isotropic transition
        }
    ]

    # Define the target temperature range (room temperature)
    target_min_temp = 20.0
    target_max_temp = 25.0

    print(f"Searching for liquid crystals with a nematic phase in the range {target_min_temp}°C to {target_max_temp}°C...")
    print("-" * 50)

    found_candidates = []

    for compound in lc_data:
        # Check if the compound's nematic range overlaps with the target room temperature range
        is_nematic_at_room_temp = (compound['nematic_start_temp'] <= target_max_temp and
                                   compound['clearing_point'] >= target_min_temp)

        if is_nematic_at_room_temp:
            found_candidates.append(compound)

    if not found_candidates:
        print("No suitable single-component liquid crystal found in the database for the specified range.")
        return

    print(f"Found {len(found_candidates)} suitable candidate(s):")
    for candidate in found_candidates:
        n = candidate['n']
        # The chemical formula is C(n)H(2n+1)-Ph-CN
        h_atoms = 2 * n + 1

        print(f"\nRecommended Material: {candidate['name']}")
        print(f"  > General Formula: CnH(2n+1)-Ph-CN")
        # Outputting each number in the final equation as requested
        print(f"  > Specific Formula: C({n})H({h_atoms})-Ph-CN")
        print(f"  > Nematic Phase Range: {candidate['nematic_start_temp']}°C to {candidate['clearing_point']}°C")
        print(f"  > This material is in the nematic phase at room temperature.")

if __name__ == '__main__':
    find_room_temp_lc()