def solve_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the maximal finite-dimensional Lie group of
    symmetries of the incompressible Navier-Stokes equations in R^3.
    """

    # The finite-dimensional group consists of the following transformations:

    # 1. Time translation: t -> t + a
    # This transformation is described by one parameter.
    dim_time_translation = 1

    # 2. Space translations: x -> x + b
    # In R^3, the translation vector b has 3 components.
    dim_space_translations = 3

    # 3. Spatial rotations: x -> Rx, u -> Ru (where R is a rotation matrix)
    # The group of rotations in 3D, SO(3), has dimension 3.
    dim_rotations = 3

    # 4. Galilean boosts: x -> x + vt, u -> u + v
    # This corresponds to changing to a frame moving with constant velocity v.
    # In R^3, the velocity vector v has 3 components.
    dim_galilean_boosts = 3

    # 5. Scaling (dilatation): (t, x, u, p) -> (c^2*t, c*x, c^-1*u, c^-2*p)
    # This transformation is described by one parameter c.
    dim_scaling = 1

    # 6. Pressure shift: p -> p + d
    # The pressure is only defined up to an additive constant.
    dim_pressure_shift = 1

    # The total dimension is the sum of the dimensions of these individual symmetry groups.
    total_dimension = (dim_time_translation +
                       dim_space_translations +
                       dim_rotations +
                       dim_galilean_boosts +
                       dim_scaling +
                       dim_pressure_shift)

    print("The dimension of the maximal finite-dimensional Lie group of symmetries of the incompressible Navier-Stokes equations in R^3 is calculated by summing the dimensions of its constituent transformations:")
    print(f"Time translation: {dim_time_translation}")
    print(f"Space translations: {dim_space_translations}")
    print(f"Rotations: {dim_rotations}")
    print(f"Galilean boosts: {dim_galilean_boosts}")
    print(f"Scaling: {dim_scaling}")
    print(f"Pressure shift: {dim_pressure_shift}")
    print(f"Total dimension = {dim_time_translation} + {dim_space_translations} + {dim_rotations} + {dim_galilean_boosts} + {dim_scaling} + {dim_pressure_shift} = {total_dimension}")

solve_navier_stokes_symmetry_dimension()