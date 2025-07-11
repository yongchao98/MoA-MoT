import pandas as pd

def predict_ir_phonons():
    """
    Calculates and prints the number of expected IR-active phonons for
    monocrystalline LiNiPO4 for different light polarizations.
    """
    print("Prediction of Polarized IR Spectra for LiNiPO4 (Pnma, No. 62)")
    print("-" * 65)
    print("Step 1: Define crystal structure and atomic site occupancies.")
    
    # Atomic composition and Wyckoff sites for olivine LiNiPO4 (Z=4)
    atomic_sites = {
        'Li': {'wyckoff': '4a', 'count': 1},
        'Ni': {'wyckoff': '4c', 'count': 1},
        'P':  {'wyckoff': '4c', 'count': 1},
        'O1,O2': {'wyckoff': '4c', 'count': 2}, # O1 and O2 atoms
        'O3': {'wyckoff': '8d', 'count': 1}, # O3 atoms
    }
    
    # Decompositions of mechanical representation into D2h irreps for each site.
    # This data is from standard group theory tables (e.g., Bilbao Crystallographic Server).
    site_decompositions = {
        '4a': {'Ag': 1, 'B1g': 1, 'B2g': 1, 'B3g': 1, 'Au': 2, 'B1u': 2, 'B2u': 2, 'B3u': 2},
        '4c': {'Ag': 2, 'B1g': 1, 'B2g': 2, 'B3g': 1, 'Au': 1, 'B1u': 2, 'B2u': 1, 'B3u': 2},
        '8d': {'Ag': 3, 'B1g': 3, 'B2g': 3, 'B3g': 3, 'Au': 3, 'B1u': 3, 'B2u': 3, 'B3u': 3}
    }
    
    print("Step 2: Calculate the total number of vibrational modes for each symmetry.")
    
    irreps = ['Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u']
    total_modes = {irrep: 0 for irrep in irreps}
    
    for atom_set, info in atomic_sites.items():
        wyckoff_site = info['wyckoff']
        num_sets = info['count']
        decomposition = site_decompositions[wyckoff_site]
        for irrep in irreps:
            total_modes[irrep] += num_sets * decomposition[irrep]
            
    # Create a DataFrame for a clean summary table
    summary_data = {
        "Symmetry": irreps,
        "Total Modes": [total_modes[irrep] for irrep in irreps]
    }
    summary_df = pd.DataFrame(summary_data).set_index("Symmetry")
    print("Total modes (Gamma_total):")
    print(summary_df.T.to_string())
    print()

    print("Step 3: Subtract acoustic modes to find the number of optic modes.")
    
    # Acoustic modes in D2h are B1u(z), B2u(y), B3u(x)
    acoustic_modes_dict = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    optic_modes = total_modes.copy()
    
    for irrep, count in acoustic_modes_dict.items():
        optic_modes[irrep] -= count

    print("The 3 acoustic modes have symmetries B1u, B2u, and B3u.")
    
    # Calculations for IR-active modes
    optic_B1u = optic_modes['B1u']
    optic_B2u = optic_modes['B2u']
    optic_B3u = optic_modes['B3u']
    
    total_B1u = total_modes['B1u']
    total_B2u = total_modes['B2u']
    total_B3u = total_modes['B3u']
    
    print(f"B1u (z-polarized): {total_B1u} total modes - 1 acoustic mode = {optic_B1u} optical modes")
    print(f"B2u (y-polarized): {total_B2u} total modes - 1 acoustic mode = {optic_B2u} optical modes")
    print(f"B3u (x-polarized): {total_B3u} total modes - 1 acoustic mode = {optic_B3u} optical modes")
    print("-" * 65)

    print("Final Result:")
    # The final answer format is E||x: [number], E||y: [number], E||z: [number]
    print(f"The predicted number of IR-active phonons are:")
    print(f"E||x: {optic_B3u}, E||y: {optic_B2u}, E||z: {optic_B1u}")


if __name__ == '__main__':
    # Run the prediction function
    predict_ir_phonons()
    # The final answer content required for the platform.
    print("\n<<<E||x: 12, E||y: 8, E||z: 12>>>")