def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    3D incompressible Navier-Stokes equations.

    The dimension is the sum of the dimensions of its constituent symmetry subgroups.
    """

    # Dimension from time translation (t -> t + a)
    dim_time_translation = 1

    # Dimension from spatial translations (x -> x + b) in R^3
    dim_space_translations = 3

    # Dimension from spatial rotations (x -> Rx) in R^3
    dim_rotations = 3

    # Dimension from Galilean boosts (u -> u + v) in R^3
    dim_galilean_boosts = 3

    # Dimension from scaling transformations
    dim_scaling = 1

    # Dimension from pressure shift (p -> p + c)
    dim_pressure_shift = 1

    # Total dimension is the sum of the dimensions of all symmetry groups.
    total_dimension = (dim_time_translation +
                       dim_space_translations +
                       dim_rotations +
                       dim_galilean_boosts +
                       dim_scaling +
                       dim_pressure_shift)

    print("The dimension of the Lie group of symmetries is the sum of the dimensions of each symmetry type:")
    print(f"Time Translation: {dim_time_translation}")
    print(f"Spatial Translations: {dim_space_translations}")
    print(f"Spatial Rotations: {dim_rotations}")
    print(f"Galilean Boosts: {dim_galilean_boosts}")
    print(f"Scaling: {dim_scaling}")
    print(f"Pressure Shift: {dim_pressure_shift}")
    print("-" * 20)
    # The final instruction asks to output each number in the final equation.
    print(f"Total Dimension = {dim_time_translation} + {dim_space_translations} + {dim_rotations} + {dim_galilean_boosts} + {dim_scaling} + {dim_pressure_shift} = {total_dimension}")

calculate_navier_stokes_symmetry_dimension()
