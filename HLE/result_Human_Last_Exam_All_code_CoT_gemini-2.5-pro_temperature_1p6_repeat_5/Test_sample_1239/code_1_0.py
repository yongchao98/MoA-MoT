def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates and prints the dimension of the Lie group of symmetries
    for the incompressible Navier-Stokes equations in R^3.
    """

    # The Lie group of symmetries for the Navier-Stokes equations is composed of
    # several types of transformations. We find the dimension of each.

    # 1. Time translation: t' = t + a. This is one parameter.
    time_translation_dim = 1

    # 2. Space translations: x' = x + v. Since space is R^3, this is a 3-vector.
    space_translation_dim = 3

    # 3. Spatial rotations: x' = R*x. The group of rotations in 3D, SO(3), has dimension 3.
    spatial_rotation_dim = 3

    # 4. Galilean boosts: x' = x - v_0*t, u' = u - v_0. The boost velocity v_0 is a 3-vector.
    galilean_boost_dim = 3

    # 5. Scaling (Dilation): t' = a^2*t, x' = a*x. This is one parameter.
    scaling_dim = 1

    # The total dimension is the sum of the dimensions of these independent generators.
    total_dimension = (
        time_translation_dim
        + space_translation_dim
        + spatial_rotation_dim
        + galilean_boost_dim
        + scaling_dim
    )

    print("The dimension is the sum of the parameters from each symmetry group:")
    print(f"- Time translation: {time_translation_dim}")
    print(f"- Space translations in R^3: {space_translation_dim}")
    print(f"- Spatial rotations in R^3: {spatial_rotation_dim}")
    print(f"- Galilean boosts in R^3: {galilean_boost_dim}")
    print(f"- Scaling (Dilation): {scaling_dim}")
    print("\nCalculating the total dimension:")

    # As requested, printing each number in the final equation
    print(f"Dimension = {time_translation_dim} + {space_translation_dim} + {spatial_rotation_dim} + {galilean_boost_dim} + {scaling_dim} = {total_dimension}")

if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()