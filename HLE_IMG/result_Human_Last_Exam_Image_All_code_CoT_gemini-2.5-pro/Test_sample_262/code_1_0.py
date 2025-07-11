def analyze_emitter_lifetimes():
    """
    Identifies which Iridium complexes are expected to have shorter lifetimes
    based on their chemical structure.
    """
    # In iridium complexes for LECs, stability is key for long lifetimes.
    # A major degradation pathway is the C-H bond ortho to the Ir-C bond.
    # Replacing this C-H with a stronger C-F bond (ortho-fluorination)
    # significantly increases stability and lifetime.

    # We represent the presence of this key feature for each complex.
    # True = has ortho-fluorination (long lifetime expected)
    # False = lacks ortho-fluorination (short lifetime expected)
    complexes = {
        1: {'has_ortho_fluorination': False},
        2: {'has_ortho_fluorination': False},
        3: {'has_ortho_fluorination': True},
        4: {'has_ortho_fluorination': True}
    }

    shorter_lifetime_complexes = []
    for num, properties in complexes.items():
        if not properties['has_ortho_fluorination']:
            shorter_lifetime_complexes.append(num)

    # Sort the list for consistent output
    shorter_lifetime_complexes.sort()

    print("Analysis based on structural stability:")
    print("Complexes with fluorination at the position ortho to the Ir-C bond are more stable and have longer lifetimes.")
    print("Complexes without this feature are less stable and have shorter lifetimes.")
    print("\nComplexes expected to show SHORTER lifetimes are:")
    
    # As requested, output each number.
    output_str = " and ".join(map(str, shorter_lifetime_complexes))
    print(f"Complex {output_str}")

analyze_emitter_lifetimes()