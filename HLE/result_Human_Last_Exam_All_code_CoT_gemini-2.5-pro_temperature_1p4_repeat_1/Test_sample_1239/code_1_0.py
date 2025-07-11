def solve_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries of the
    3D incompressible Navier-Stokes equations.
    """

    # The symmetry group of the 3D Navier-Stokes equations corresponds to the
    # Schrödinger group sch(3). We sum the number of generators for each
    # independent symmetry transformation.

    # Number of generators for time translation
    time_translations = 1

    # Number of generators for space translations in R^3
    space_translations = 3

    # Number of generators for space rotations (SO(3))
    space_rotations = 3

    # Number of generators for Galilean boosts in R^3
    galilean_boosts = 3

    # Number of generators for dilations (scaling)
    dilations = 1

    # Number of generators for special conformal transformations
    special_conformal_transformations = 1

    # The total dimension is the sum of the individual generators.
    total_dimension = (time_translations +
                       space_translations +
                       space_rotations +
                       galilean_boosts +
                       dilations +
                       special_conformal_transformations)

    print("To find the dimension of the Lie group of symmetries, we sum the number of generators for each type of transformation:")
    print(f"1. Time translations: {time_translations}")
    print(f"2. Space translations: {space_translations}")
    print(f"3. Space rotations: {space_rotations}")
    print(f"4. Galilean boosts: {galilean_boosts}")
    print(f"5. Dilations (scaling): {dilations}")
    print(f"6. Special conformal transformations: {special_conformal_transformations}")
    print("-" * 30)
    
    # Final equation showing the sum
    print("Total Dimension =",
          f" {time_translations} + {space_translations} + {space_rotations} + "
          f"{galilean_boosts} + {dilations} + {special_conformal_transformations} =",
          f" {total_dimension}")

solve_navier_stokes_symmetry_dimension()