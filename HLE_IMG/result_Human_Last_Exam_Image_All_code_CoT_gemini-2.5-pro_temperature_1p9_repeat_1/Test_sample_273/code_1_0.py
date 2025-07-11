def determine_molecule_configuration():
    """
    This script outlines the reasoning for determining the absolute configuration
    of the chiral centers in the provided molecule.
    """
    print("Step-by-step determination of absolute configuration:")
    print("IUPAC Name: 4-amino-3-ethyl-4-methylpentan-2-ol")
    print("Chiral centers: C2, C3, C4")
    print("-" * 40)

    # --- C3 Analysis ---
    c3_config = 'R'
    print("Analysis of C3:")
    print("  Priorities: 1:-CH(OH)CH3 > 2:-CH(...) > 3:-CH2CH3 > 4:-H")
    print("  Geometry: Lowest priority group (-H) is pointing away (dash).")
    print("  Sequence 1 -> 2 -> 3 is clockwise.")
    print(f"  Configuration of C3 is {c3_config}.")
    print("-" * 40)

    # --- C4 Analysis ---
    c4_config = 'R'
    print("Analysis of C4:")
    print("  Priorities: 1:-CH2NH2 > 2:-C(...) > 3:-CH3 > 4:-H")
    print("  Geometry: Lowest priority group (-H) is pointing forward (wedge).")
    print("  Sequence 1 -> 2 -> 3 is counter-clockwise (S).")
    print("  Result is reversed because -H is forward.")
    print(f"  Configuration of C4 is {c4_config}.")
    print("-" * 40)

    # --- C2 Analysis ---
    c2_config = 'R'
    print("Analysis of C2:")
    print("  Priorities: 1:-OH > 2:-C(...) > 3:-CH3 > 4:-H")
    print("  Geometry: Lowest priority group (-H) is pointing forward (wedge).")
    print("  Sequence 1 -> 2 -> 3 is counter-clockwise (S).")
    print("  Result is reversed because -H is forward.")
    print(f"  Configuration of C2 is {c2_config}.")
    print("-" * 40)

    # --- Final Result ---
    print("Final Absolute Configuration:")
    print(f"C2 is {c2_config}")
    print(f"C3 is {c3_config}")
    print(f"C4 is {c4_config}")
    print("\nSummary:")
    # Printing each part of the final configuration string
    print(f"The first chiral center is at carbon 2 and has configuration {c2_config}")
    print(f"The second chiral center is at carbon 3 and has configuration {c3_config}")
    print(f"The third chiral center is at carbon 4 and has configuration {c4_config}")


determine_molecule_configuration()