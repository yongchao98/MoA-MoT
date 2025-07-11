def solve_navier_stokes_symmetry_dimension():
    """
    Calculates and prints the dimension of the Lie group of symmetries
    for the 3D incompressible Navier-Stokes equations.
    """
    
    # The dimension of the Lie group is the number of independent parameters
    # in the continuous transformations that leave the equations invariant.
    # We sum the dimensions from each type of symmetry.

    # 1. Time translation: t -> t + s (1 parameter: s)
    time_translation = 1

    # 2. Space translations: x -> x + a (3 parameters: a_x, a_y, a_z)
    space_translations = 3

    # 3. Space rotations: x -> Rx (3 parameters: e.g., Euler angles)
    space_rotations = 3

    # 4. Galilean boosts: x -> x + ut, v -> v + u (3 parameters: u_x, u_y, u_z)
    galilean_boosts = 3

    # 5. Scaling (dilatation): t -> λ^2 t, x -> λx, etc. (1 parameter: λ)
    scaling = 1

    # 6. Pressure shift: p -> p + c (1 parameter: c)
    pressure_shift = 1

    # The total dimension is the sum of these individual dimensions.
    total_dimension = (
        time_translation +
        space_translations +
        space_rotations +
        galilean_boosts +
        scaling +
        pressure_shift
    )

    print("The total dimension of the Lie group of symmetries is the sum of dimensions from each type of symmetry:")
    print(f"Dimension = "
          f"{time_translation} (Time Translation) + "
          f"{space_translations} (Space Translations) + "
          f"{space_rotations} (Space Rotations) + "
          f"{galilean_boosts} (Galilean Boosts) + "
          f"{scaling} (Scaling) + "
          f"{pressure_shift} (Pressure Shift)")
    
    print(f"Total Dimension = {total_dimension}")

if __name__ == "__main__":
    solve_navier_stokes_symmetry_dimension()