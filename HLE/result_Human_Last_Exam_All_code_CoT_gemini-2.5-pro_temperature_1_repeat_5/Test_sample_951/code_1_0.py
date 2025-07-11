import numpy as np

def solve_phonon_modes():
    """
    Calculates the number of IR-active phonon modes for LiNiPO4 (Pnma space group)
    using group theory (Factor Group Analysis).
    """

    # 1. Define the character table for the D2h point group.
    # Operations: E, C2(z), C2(y), C2(x), i, sigma(xy), sigma(xz), sigma(yz)
    D2h_char_table = {
        'Ag':  [1, 1, 1, 1, 1, 1, 1, 1],
        'B1g': [1, 1, -1, -1, 1, 1, -1, -1],
        'B2g': [1, -1, 1, -1, 1, -1, 1, -1],
        'B3g': [1, -1, -1, 1, 1, -1, -1, 1],
        'Au':  [1, 1, 1, 1, -1, -1, -1, -1],
        'B1u': [1, 1, -1, -1, -1, -1, 1, 1],  # Transforms as z -> IR active E||z
        'B2u': [1, -1, 1, -1, -1, 1, -1, 1],  # Transforms as y -> IR active E||y
        'B3u': [1, -1, -1, 1, -1, 1, 1, -1],  # Transforms as x -> IR active E||x
    }
    
    # Order of the point group
    group_order = 8.0

    # 2. Information for LiNiPO4 in space group Pnma (No. 62)
    # Total atoms in primitive cell = 4 * (Li+Ni+P+4O) = 28 atoms.
    # The number of atoms left unmoved by the symmetry operations (R) of D2h.
    # R = [E, C2(z), C2(y), C2(x), i, sigma(xy), sigma(xz), sigma(yz)]
    # For Pnma: Li(4a), Ni(4c), P(4c), O1(4c), O2(4c), O3(8d)
    # n(E)=28, n(i)=4 (Li atoms), n(sigma(xz))=16 (Ni,P,O1,O2 atoms on mirror plane)
    # All other operations leave 0 atoms unmoved.
    unmoved_atoms = np.array([28, 0, 0, 0, 4, 0, 16, 0])

    # 3. Characters of the vector representation for a general vector (x,y,z)
    # chi_v(R) = 1 + 2*cos(theta) for pure rotations
    # chi_v(R) = -1 + 2*cos(theta) for improper rotations (roto-inversions)
    # R = [E, C2, C2, C2, i, sigma, sigma, sigma]
    # E(theta=0), C2(theta=180), i(S2, theta=180), sigma(S1, theta=0)
    char_vector = np.array([3, -1, -1, -1, -3, 1, 1, 1])

    # 4. Calculate the characters of the total reducible representation
    # chi_total(R) = n(R) * chi_vector(R)
    chi_total = unmoved_atoms * char_vector

    # 5. Decompose the total representation into irreducible representations
    optical_modes = {}
    for irrep, characters in D2h_char_table.items():
        n_irrep = (1 / group_order) * np.sum(np.array(characters) * chi_total)
        optical_modes[irrep] = int(round(n_irrep))

    # 6. Subtract acoustic modes to get optical modes
    # The acoustic modes transform as x, y, z, which are B3u, B2u, and B1u.
    acoustic_irreps = ['B1u', 'B2u', 'B3u']
    for irrep in acoustic_irreps:
        if irrep in optical_modes:
            optical_modes[irrep] -= 1
            
    # 7. Extract the numbers for IR-active modes and print the result
    num_x = optical_modes.get('B3u', 0)
    num_y = optical_modes.get('B2u', 0)
    num_z = optical_modes.get('B1u', 0)

    # The problem asks to output each number in the final equation.
    # The final equation is the formatted string below.
    print(f"E||x: {num_x}, E||y: {num_y}, E||z: {num_z}")

solve_phonon_modes()