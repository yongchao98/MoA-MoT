def analyze_crystal_classes():
    """
    Analyzes crystal classes based on optical properties to answer the user's question.
    """
    # Data for crystal classes. A class is optically active if and only if it's chiral.
    # Hermann-Mauguin symbols are used. '4m' is an alternative for '4mm'.
    crystal_data = {
        # Achiral & Polar (can show directional optical activity)
        'm':    {'chiral': False, 'polar': True},
        'mm2':  {'chiral': False, 'polar': True},
        '3m':   {'chiral': False, 'polar': True},
        '4mm':  {'chiral': False, 'polar': True},
        '6mm':  {'chiral': False, 'polar': True},
        # Achiral & Non-polar
        '-1':   {'chiral': False, 'polar': False},
        '2/m':  {'chiral': False, 'polar': False},
        'mmm':  {'chiral': False, 'polar': False},
        '-4':   {'chiral': False, 'polar': False},
        '4/m':  {'chiral': False, 'polar': False},
        '-42m': {'chiral': False, 'polar': False},
        '-3':   {'chiral': False, 'polar': False},
        '-3m':  {'chiral': False, 'polar': False},
        '-6':   {'chiral': False, 'polar': False},
        '-62m': {'chiral': False, 'polar': False}, # same as -6m2
        '-43m': {'chiral': False, 'polar': False},
        # Chiral (and therefore Optically Active)
        '1':    {'chiral': True, 'polar': True},
        '2':    {'chiral': True, 'polar': True},
        '3':    {'chiral': True, 'polar': True},
        '4':    {'chiral': True, 'polar': True},
        '6':    {'chiral': True, 'polar': True},
    }

    print("--- Step 1: Analysis of the Question's Premise ---")
    print("The question asks for crystal classes with three properties:")
    print("1. Achiral (possesses a mirror plane or inversion center)")
    print("2. Non-polar (has no net electric dipole moment)")
    print("3. Optically Active (rotates the plane of polarized light)")
    print("\nPrinciple of Physics: Natural optical activity is only possible in chiral crystal classes.")
    print("Conclusion: The conditions 'achiral' and 'optically active' are mutually exclusive. Therefore, no crystal class satisfies the question as written.")

    print("\n--- Step 2: Searching for Classes that fit the Impossible Criteria ---")
    impossible_classes = []
    for name, props in crystal_data.items():
        # is_optically_active is equivalent to props['chiral']
        if not props['chiral'] and not props['polar'] and props['chiral']:
            impossible_classes.append(name)
    
    # The final equation demonstrates the contradiction
    print("Final Equation (original query):")
    print("Find(Class) WHERE IsAchiral = True AND IsNonPolar = True AND IsOpticallyActive = True")
    # This is the key part of the 'equation' showing the logic
    print("Let IsOpticallyActive = IsChiral. The equation becomes:")
    print("Find(Class) WHERE IsChiral = False AND IsNonPolar = True AND IsChiral = True")
    print(f"Result: {impossible_classes} (The set is empty, as expected).")

    print("\n--- Step 3: Proposing a Correction to find the Intended Answer ---")
    print("A plausible interpretation is that the question contains a typo, and 'non-polar' should be 'polar'.")
    print("The corrected question would be: 'What achiral and POLAR crystal classes...'")
    print("This might refer to directional optical activity, a phenomenon possible in these classes.")
    
    corrected_classes = []
    for name, props in crystal_data.items():
        if not props['chiral'] and props['polar']:
            corrected_classes.append(name)

    print("\nFinal Equation (corrected query):")
    print("Find(Class) WHERE IsAchiral = True AND IsPolar = True")
    print(f"Result: {sorted(corrected_classes)}")
    
    print("\n--- Step 4: Analyzing Answer Choices ---")
    print("We now check which answer choice matches the result from the corrected query.")
    choices = {
        'A': ['m', 'mm2'],
        'B': ['-6', '-62m', '-43m'],
        'C': ['3m', '4m', '6mm'],
        'D': ['-4', '-42m'],
        'E': ['1', '2', '3', '4', '6']
    }
    print(f"Choice C, which is {choices['C']} (treating '4m' as '4mm'), is a perfect subset of the result from our corrected query.")

analyze_crystal_classes()