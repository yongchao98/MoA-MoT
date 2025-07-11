def find_unstable_complexes():
    """
    Identifies iridium complexes with structural features known to cause short lifetimes in LECs.

    The key principle is that a fluorine atom (F) on the phenylpyridine ligand at a position
    'ortho' (adjacent) to the Iridium-Carbon (Ir-C) bond leads to molecular instability
    and thus a shorter device lifetime.
    """

    # We represent our analysis of each complex in a dictionary.
    # The value 'True' means the complex has an unstable ortho-C-F bond.
    complex_analysis = {
        'Complex 1': True,
        'Complex 2': False,
        'Complex 3': True,
        'Complex 4': True
    }

    short_lifetime_complexes = []
    print("Analysis of Complex Stability:")
    print("-" * 30)

    for complex_name, has_ortho_F in complex_analysis.items():
        complex_number = int(complex_name.split()[-1])
        if has_ortho_F:
            stability = "Unstable (shorter lifetime expected)"
            short_lifetime_complexes.append(complex_number)
        else:
            stability = "Stable (longer lifetime expected)"
        
        print(f"{complex_name}: Has ortho-C-F bond? {has_ortho_F}. --> {stability}")

    short_lifetime_complexes.sort()
    
    # Format the final list of numbers as required
    output_numbers = ", ".join(map(str, short_lifetime_complexes))

    print("\nConclusion:")
    print("The complexes expected to show shorter lifetimes are those with an ortho-C-F bond.")
    print(f"These are complexes: {output_numbers}")

# Execute the function to get the result
find_unstable_complexes()