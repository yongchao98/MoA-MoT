def find_unstable_complexes():
    """
    Analyzes the Iridium complexes to predict which ones have shorter lifetimes in LECs.

    The stability of these Ir(III) complexes is a key factor for device lifetime.
    A major degradation pathway is the dissociation of the cyclometalating phenylpyridine ligand.
    The strength of the Iridium-Carbon (Ir-C) bond is crucial.

    Structural Analysis:
    - All complexes have the same ancillary bipyridine ligand.
    - The difference lies in the fluorine substituents on the phenylpyridine ligands.
    - Electronic Effect: Fluorine is electron-withdrawing, which generally strengthens the Ir-C bond and increases stability.
    - Steric Effect: Substituents in the ortho position (2'-position) to the Ir-C bond cause steric hindrance. This strain weakens the Ir-C bond, leading to lower stability and faster degradation.

    - Complex 1: Fluorines at 4', 5' positions. No steric hindrance. Expected to be stable.
    - Complex 2: Fluorine at 4' position. No steric hindrance. Expected to be stable.
    - Complex 3: Fluorine at 2' (ortho) position. Significant steric hindrance. Expected to be unstable.
    - Complex 4: Fluorines at 2' (ortho) and 4' positions. Significant steric hindrance from the 2'-F. Expected to be unstable.

    Conclusion: Complexes with ortho-fluorine substituents (Complexes 3 and 4) are expected to be less stable and have shorter lifetimes.
    """
    complexes = [
        {'id': 1, 'ortho_F': False, 'comment': 'F at 4\', 5\'. Stable.'},
        {'id': 2, 'ortho_F': False, 'comment': 'F at 4\'. Stable.'},
        {'id': 3, 'ortho_F': True, 'comment': 'F at 2\'. Steric hindrance. Unstable.'},
        {'id': 4, 'ortho_F': True, 'comment': 'F at 2\', 4\'. Steric hindrance. Unstable.'}
    ]

    shorter_lifetime_complexes = []
    for c in complexes:
        if c['ortho_F']:
            shorter_lifetime_complexes.append(c['id'])

    # Format the output string
    complex_list_str = ', '.join(map(str, shorter_lifetime_complexes))
    print(f"Complexes with fluorine atoms at the ortho-position (2') of the phenyl ring experience significant steric hindrance.")
    print("This steric strain weakens the Ir-C bond, leading to lower stability and shorter device lifetimes.")
    print(f"Therefore, complexes [{complex_list_str}] are expected to show shorter lifetimes.")

find_unstable_complexes()