def calculate_helical_pattern():
    """
    Calculates the predicted H-bond ring size for an alternating
    alanine/epsilon-amino acid foldamer.
    """
    # Step 1: Define monomer properties
    # n_full: total number of backbone atoms (N, Ca, ..., C')
    # n_part: number of backbone atoms from Ca to C'
    monomers = {
        'alanine': {'name': 'Alanine (alpha-AA)', 'n_full': 3, 'n_part': 2},
        'epsilon': {'name': 'Epsilon-AA', 'n_full': 7, 'n_part': 6}
    }

    print("--- Monomer Properties ---")
    print(f"Alanine backbone atoms: {monomers['alanine']['n_full']}")
    print(f"Epsilon-AA backbone atoms: {monomers['epsilon']['n_full']}\n")

    # Step 2 & 3: Calculate ring size for the i -> i+2 H-bond pattern
    print("--- Calculating H-bond Ring Size for i -> i+2 Pattern ---")
    
    # Case 1: Ala -> Ala bond, spanning an Epsilon residue
    res_i_1 = 'alanine'
    spanned_1 = ['epsilon']
    # Formula: Ring Size = 2 + N_part(i) + N_full(spanned i+1)
    ring_size_1 = 2 + monomers[res_i_1]['n_part'] + monomers[spanned_1[0]]['n_full']
    
    print("Case 1: Ala(i) -> Ala(i+2) bond spanning an Epsilon-AA residue")
    print(f"Ring Size = 2 + {monomers[res_i_1]['n_part']} (from Ala) + {monomers[spanned_1[0]]['n_full']} (from Eps)")
    print(f"Resulting Ring Size = {ring_size_1}\n")

    # Case 2: Eps -> Eps bond, spanning an Alanine residue
    res_i_2 = 'epsilon'
    spanned_2 = ['alanine']
    # Formula: Ring Size = 2 + N_part(i) + N_full(spanned i+1)
    ring_size_2 = 2 + monomers[res_i_2]['n_part'] + monomers[spanned_2[0]]['n_full']

    print("Case 2: Eps(i) -> Eps(i+2) bond spanning an Alanine residue")
    print(f"Ring Size = 2 + {monomers[res_i_2]['n_part']} (from Eps) + {monomers[spanned_2[0]]['n_full']} (from Ala)")
    print(f"Resulting Ring Size = {ring_size_2}\n")
    
    # Step 4 & 5: Conclusion
    print("--- Conclusion ---")
    print(f"Both calculations for the i->i+2 H-bond pattern predict a ring of size {ring_size_1}.")
    print("This suggests the fundamental structure is an 11-helix.")
    print("The notation x/y suggests alternating ring sizes, likely a perturbation of the base 11-helix.")
    print("Of the choices containing 11 (11/9 and 11/13), the 11/13 pattern is most plausible as it combines two known stable H-bond motifs.")

calculate_helical_pattern()