def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the maximal finite-dimensional Lie group of
    symmetries of the incompressible Navier-Stokes equations in R^3.

    The total dimension is the sum of the dimensions of the individual known
    symmetry groups.
    """

    # 1. Time translation: Shifting the time origin (t -> t + a).
    # This is a one-parameter group, so its dimension is 1.
    time_translation_dim = 1

    # 2. Spatial translations: Shifting the spatial origin (x -> x + b).
    # In R^3, the translation vector b has 3 components. Dimension is 3.
    space_translation_dim = 3

    # 3. Spatial rotations: Rotating the coordinate system (x -> R * x).
    # The group of rotations in R^3, SO(3), has dimension 3.
    rotation_dim = 3

    # 4. Galilean boosts: Transforming to a constantly moving frame.
    # The boost velocity vector v has 3 components. Dimension is 3.
    galilean_boost_dim = 3

    # 5. Scaling transformations: Rescaling space, time, velocity, and pressure.
    # This is a one-parameter group (t -> a^2*t, x -> a*x, ...). Dimension is 1.
    scaling_dim = 1

    # The total dimension is the sum of these individual dimensions.
    total_dimension = (
        time_translation_dim
        + space_translation_dim
        + rotation_dim
        + galilean_boost_dim
        + scaling_dim
    )

    # Print the breakdown of the calculation as an equation.
    print("The total dimension is the sum of the dimensions of each symmetry component:")
    print(
        f"Total Dimension = {time_translation_dim} (time translation) "
        f"+ {space_translation_dim} (space translations) "
        f"+ {rotation_dim} (rotations) "
        f"+ {galilean_boost_dim} (Galilean boosts) "
        f"+ {scaling_dim} (scaling)"
    )
    print(f"Total Dimension = {total_dimension}")

if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()