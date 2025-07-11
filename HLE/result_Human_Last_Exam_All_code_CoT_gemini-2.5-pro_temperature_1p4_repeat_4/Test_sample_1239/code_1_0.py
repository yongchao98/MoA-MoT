def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries of the
    incompressible Navier-Stokes equations in R^3.

    The dimension is found by summing the number of infinitesimal generators
    corresponding to each family of symmetries.
    """

    # 1. Time translations: t' = t + s
    # The equations are autonomous.
    time_translations = 1

    # 2. Space translations: x' = x + a
    # The equations are spatially homogeneous. R^3 gives 3 generators.
    space_translations = 3

    # 3. Spatial rotations: x' = R * x
    # The equations are isotropic. The rotation group SO(3) is 3-dimensional.
    spatial_rotations = 3

    # 4. Galilean boosts: x' = x + v*t, u' = u + v
    # Invariance under change of inertial frame. 3 generators for velocity v in R^3.
    galilean_boosts = 3

    # 5. Scaling (Dilation) symmetry:
    # A specific scaling of time, space, velocity, and pressure leaves the equations invariant.
    # This provides a single generator.
    scaling_symmetries = 1

    # 6. Pressure shift: p' = p + c
    # The equations only depend on the gradient of pressure, so a constant shift is a symmetry.
    pressure_shift = 1

    # The dimension of the Lie algebra (and thus the Lie group) is the sum
    # of the number of generators from each independent family of symmetries.
    total_dimension = (time_translations +
                       space_translations +
                       spatial_rotations +
                       galilean_boosts +
                       scaling_symmetries +
                       pressure_shift)

    print("The dimension of the Lie group of symmetries of the incompressible Navier-Stokes equations in R^3 is the sum of the number of generators from each symmetry family:")
    print(f"  - Time translations:      {time_translations}")
    print(f"  - Space translations:     {space_translations}")
    print(f"  - Spatial rotations:      {spatial_rotations}")
    print(f"  - Galilean boosts:        {galilean_boosts}")
    print(f"  - Scaling symmetry:       {scaling_symmetries}")
    print(f"  - Constant pressure shift:{pressure_shift}")
    print("-" * 50)
    print(f"Total dimension = {time_translations} + {space_translations} + {spatial_rotations} + {galilean_boosts} + {scaling_symmetries} + {pressure_shift} = {total_dimension}")

calculate_navier_stokes_symmetry_dimension()