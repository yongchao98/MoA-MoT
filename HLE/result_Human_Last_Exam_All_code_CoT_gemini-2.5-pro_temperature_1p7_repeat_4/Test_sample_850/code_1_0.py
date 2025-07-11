def solve_crystal_class_query():
    """
    Analyzes the 32 crystal classes to find which are achiral, non-polar,
    and optically active, explaining the relationships between these properties.
    """

    # The 11 chiral crystal classes. These are the ONLY classes that are optically active.
    # Optical activity is a direct consequence of chirality.
    optically_active_classes = sorted(['1', '2', '222', '3', '32', '4', '422', '6', '622', '23', '432'])

    # The 10 polar crystal classes. These exhibit properties like pyroelectricity.
    polar_classes = sorted(['1', '2', '3', '4', '6', 'm', 'mm2', '3m', '4mm', '6mm'])

    # All 32 crystallographic point groups (crystal classes)
    all_classes = sorted([
        '1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m',
        '422', '4mm', '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6',
        '-6', '6/m', '622', '6mm', '-6m2', '6/mmm', '23', 'm-3', '432',
        '-43m', 'm-3m'
    ])

    # Find the achiral and non-polar classes based on the definitions.
    achiral_non_polar_classes = []
    for c in all_classes:
        is_achiral = c not in optically_active_classes
        is_non_polar = c not in polar_classes
        if is_achiral and is_non_polar:
            achiral_non_polar_classes.append(c)

    # Attempt to find classes that meet the user's contradictory criteria.
    # is_optically_active == is_chiral, so (is_achiral and is_optically_active) is impossible.
    contradictory_classes = []
    for c in all_classes:
        is_achiral = c not in optically_active_classes
        is_non_polar = c not in polar_classes
        is_optically_active = c in optically_active_classes # This is the key definition
        if is_achiral and is_non_polar and is_optically_active:
            contradictory_classes.append(c)
    
    # --- Output the Explanation and Results ---
    print("--- Analysis of Crystal Class Properties ---")
    print("\n[Principle]")
    print("A crystal class has the correct symmetry for optical activity if and only if it is chiral.")
    print("A chiral class lacks a center of inversion and a mirror plane.")
    print("The term 'achiral' means 'not chiral'.")
    print("\n[Conclusion]")
    print("Therefore, it is impossible for a crystal class to be both 'achiral' and 'optically active'.")
    print("Your query asks for classes that are simultaneously achiral and optically active, which is a contradiction.")
    
    print("\n--- Search Results ---")
    if not contradictory_classes:
        print("Achiral, non-polar, and optically active classes found: None")
    else:
        # This part of the code will never be reached, but is included for logical completeness.
        print("Achiral, non-polar, and optically active classes found: ", ', '.join(contradictory_classes))

    print("\n--- For Reference ---")
    print("Here are the correct lists for the properties you mentioned:")
    
    print("\n1. Optically Active (i.e., Chiral) Crystal Classes (11 total):")
    print(', '.join(optically_active_classes))

    print("\n2. Achiral AND Non-Polar Crystal Classes (15 total):")
    print(', '.join(achiral_non_polar_classes))


if __name__ == '__main__':
    solve_crystal_class_query()