def solve_nanoribbon_classification():
    """
    Analyzes and classifies three graphene nanoribbon band structures.
    """

    # --- Ribbon 1 Analysis ---
    # Edge: Armchair (A) - symmetric bands, no flat edge states at E=0.
    # Width: 8 - counted 8 bands above the Fermi level.
    # Band: Semiconducting (1) - has a clear band gap at E=0.
    edge_1 = 'A'
    width_1 = 8
    band_1 = 1
    class_1 = f"{edge_1}{width_1}{band_1}"

    # --- Ribbon 2 Analysis ---
    # Edge: Zigzag (Z) - characteristic flat bands at E=0.
    # Width: 7 - counted 7 bands above the Fermi level.
    # Band: Metallic (0) - bands cross the Fermi level.
    edge_2 = 'Z'
    width_2 = 7
    band_2 = 0
    class_2 = f"{edge_2}{width_2}{band_2}"

    # --- Ribbon 3 Analysis ---
    # Edge: Armchair (A) - symmetric bands, no flat edge states.
    # Width: 5 - counted 5 bands above the Fermi level.
    # Band: Metallic (0) - bands touch at E=0 (zero-gap).
    edge_3 = 'A'
    width_3 = 5
    band_3 = 0
    class_3 = f"{edge_3}{width_3}{band_3}"

    # Concatenate the classifications for all three ribbons
    final_classification = class_1 + class_2 + class_3
    
    # The prompt asks to "output each number in the final equation".
    # We interpret this as showing how the final string is constructed.
    # The classification for Ribbon 1 is Edge='A', Width=8, Band=1.
    # The classification for Ribbon 2 is Edge='Z', Width=7, Band=0.
    # The classification for Ribbon 3 is Edge='A', Width=5, Band=0.
    # Concatenating them gives: A81Z70A50
    print(final_classification)

solve_nanoribbon_classification()