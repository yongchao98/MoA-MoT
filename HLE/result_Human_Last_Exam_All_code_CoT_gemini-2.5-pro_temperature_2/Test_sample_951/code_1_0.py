import collections

def calculate_ir_phonons():
    """
    Performs a factor group analysis for LiNiPO4 (space group Pnma) to determine
    the number of IR-active phonons for each polarization.
    """

    # Irreducible representations for the D2h point group
    irreps = ['Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u']

    # Decomposition of vibrational modes for atoms at Wyckoff positions 4c and 8d in Pnma
    # Source: Bilbao Crystallographic Server
    gamma_4c = {'Ag': 1, 'B1g': 1, 'B2g': 1, 'B3g': 2, 'Au': 2, 'B1u': 1, 'B2u': 2, 'B3u': 1}
    gamma_8d = {'Ag': 3, 'B1g': 3, 'B2g': 3, 'B3g': 3, 'Au': 3, 'B1u': 3, 'B2u': 3, 'B3u': 3}

    # Number of atom sets on each Wyckoff site for LiNiPO4
    # 5 sets on 4c: {Li}, {Ni}, {P}, {O1}, {O2}
    # 1 set on 8d: {O3}
    num_4c_sets = 5
    num_8d_sets = 1

    print("--- Step 1: Calculate Total Mechanical Representation (Gamma_total) ---")
    print(f"Gamma_total = {num_4c_sets} * Gamma_4c + {num_8d_sets} * Gamma_8d\n")

    gamma_total = collections.defaultdict(int)
    for ir in irreps:
        contribution_4c = num_4c_sets * gamma_4c.get(ir, 0)
        contribution_8d = num_8d_sets * gamma_8d.get(ir, 0)
        gamma_total[ir] = contribution_4c + contribution_8d
        print(f"  {ir}: {num_4c_sets} * {gamma_4c.get(ir, 0)} + {num_8d_sets} * {gamma_8d.get(ir, 0)} = {gamma_total[ir]}")
        
    print("\nResulting Gamma_total:")
    total_modes_str = " + ".join([f"{count}{ir}" for ir, count in sorted(gamma_total.items())])
    print(f"Gamma_total = {total_modes_str}\n")
    
    # Subtract acoustic modes
    print("--- Step 2: Subtract Acoustic Modes to find Optical Modes ---")
    gamma_optical = gamma_total.copy()
    acoustic_modes = {'B3u': 1, 'B2u': 1, 'B1u': 1} # Tx, Ty, Tz respectively

    print("Gamma_acoustic = B1u + B2u + B3u\n")
    print("Gamma_optical = Gamma_total - Gamma_acoustic\n")
    for ir, count in acoustic_modes.items():
        print(f"  Subtracting {count} {ir} mode:")
        print(f"    {ir}: {gamma_optical[ir]} - {count} = {gamma_optical[ir] - count}")
        gamma_optical[ir] -= count

    print("\nResulting Gamma_optical:")
    optical_modes_str = " + ".join([f"{count}{ir}" for ir, count in sorted(gamma_optical.items())])
    print(f"Gamma_optical = {optical_modes_str}\n")
    
    # Extract IR active modes
    # B3u is active for E||x
    # B2u is active for E||y
    # B1u is active for E||z
    num_x = gamma_optical['B3u']
    num_y = gamma_optical['B2u']
    num_z = gamma_optical['B1u']
    
    print("--- Step 3: Identify IR Active Phonons for each Polarization ---")
    print("IR activity corresponds to B1u(z), B2u(y), and B3u(x) modes.")
    print("\nFinal Count:")
    print(f"E||x: {num_x} ({num_x} B3u modes)")
    print(f"E||y: {num_y} ({num_y} B2u modes)")
    print(f"E||z: {num_z} ({num_z} B1u modes)")
    
    print("\nFinal Answer Formatted:")
    print(f"E||x: [{num_x}], E||y: [{num_y}], E||z: [{num_z}]")

calculate_ir_phonons()