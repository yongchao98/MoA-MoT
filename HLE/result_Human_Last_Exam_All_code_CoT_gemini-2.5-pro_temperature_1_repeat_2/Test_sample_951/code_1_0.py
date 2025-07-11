import collections

def solve_phonon_modes():
    """
    Calculates the number of IR-active phonon modes for LiNiPO4.
    """
    # D2h character table and basis vectors
    d2h_table = {
        'Ag':  [1, 1, 1, 1, 1, 1, 1, 1],
        'B1g': [1, 1, -1, -1, 1, 1, -1, -1],
        'B2g': [1, -1, 1, -1, 1, -1, 1, -1],
        'B3g': [1, -1, -1, 1, 1, -1, -1, 1],
        'Au':  [1, 1, 1, 1, -1, -1, -1, -1],
        'B1u': [1, 1, -1, -1, -1, -1, 1, 1],
        'B2u': [1, -1, 1, -1, -1, 1, -1, 1],
        'B3u': [1, -1, -1, 1, -1, 1, 1, -1],
    }
    # Symmetry operations for D2h point group
    # Order: E, C2(z), C2(y), C2(x), i, sigma(xy), sigma(xz), sigma(yz)
    operations = ['E', 'C2(z)', 'C2(y)', 'C2(x)', 'i', 'sigma(xy)', 'sigma(xz)', 'sigma(yz)']
    group_order = len(operations)

    # Characters of the vector representation for each operation
    char_vector = [3, -1, -1, -1, -3, 1, 1, 1]

    # Number of atoms unshifted by each operation for LiNiPO4 in Pnma
    # E: All 28 atoms
    # i: 4 Li atoms on inversion centers
    # sigma(xz): 16 atoms (Ni, P, 2xO) on mirror planes
    # Other ops correspond to screw/glide planes, so 0 atoms are unshifted.
    n_unshifted = [28, 0, 0, 0, 4, 0, 16, 0]

    print("Step 1: Calculating characters for the total representation (Gamma_total)")
    # Character of the total representation: chi_total = n_unshifted * chi_vector
    chi_total = [n * v for n, v in zip(n_unshifted, char_vector)]
    print(f"Characters: {chi_total}\n")


    print("Step 2: Decomposing Gamma_total into irreducible representations")
    total_modes = collections.OrderedDict()
    for irrep, characters in d2h_table.items():
        count = sum(c * ct for c, ct in zip(characters, chi_total)) // group_order
        if count > 0:
            total_modes[irrep] = count

    # Print the decomposition of total modes
    gamma_total_str = " + ".join([f"{count}{irrep}" for irrep, count in total_modes.items()])
    print(f"Gamma_total = {gamma_total_str}\n")


    print("Step 3: Subtracting acoustic modes to find optical modes (Gamma_optic)")
    # Acoustic modes in D2h transform as B1u(z), B2u(y), B3u(x)
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    print("Acoustic modes are: 1B1u + 1B2u + 1B3u\n")
    
    optical_modes = total_modes.copy()
    for irrep, count in acoustic_modes.items():
        optical_modes[irrep] -= count

    # Print the decomposition of optical modes
    gamma_optic_str = " + ".join([f"{count}{irrep}" for irrep, count in optical_modes.items()])
    print(f"Gamma_optic = {gamma_optic_str}\n")


    print("Step 4: Identifying IR-active modes by polarization")
    # IR active modes correspond to B1u, B2u, B3u
    ir_active_modes = {
        'x': optical_modes.get('B3u', 0), # B3u transforms as x
        'y': optical_modes.get('B2u', 0), # B2u transforms as y
        'z': optical_modes.get('B1u', 0)  # B1u transforms as z
    }

    print("The number of expected phonons in each polarization is:")
    print(f"E||x: {ir_active_modes['x']}, E||y: {ir_active_modes['y']}, E||z: {ir_active_modes['z']}")

solve_phonon_modes()