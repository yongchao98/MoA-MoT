def calculate_lie_group_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    3D incompressible Navier-Stokes equations.
    """

    # Dimension from time translation (t -> t + a)
    time_translation_dim = 1

    # Dimension from space translations (x -> x + b) in R^3
    space_translation_dim = 3

    # Dimension from spatial rotations (x -> Rx) in R^3
    rotation_dim = 3

    # Dimension from Galilean boosts (x -> x - vt, u -> u - v) in R^3
    galilean_boost_dim = 3

    # Dimension from scaling symmetry (t -> a^2*t, x -> a*x, ...)
    scaling_dim = 1

    # Total dimension is the sum of all individual symmetry dimensions
    total_dimension = (
        time_translation_dim
        + space_translation_dim
        + rotation_dim
        + galilean_boost_dim
        + scaling_dim
    )

    print("The dimension of the Lie group of symmetries for the 3D incompressible Navier-Stokes equations is calculated by summing the dimensions of its constituent symmetries:")
    print(f"Time Translations: {time_translation_dim}")
    print(f"Space Translations: {space_translation_dim}")
    print(f"Rotations: {rotation_dim}")
    print(f"Galilean Boosts: {galilean_boost_dim}")
    print(f"Scaling Symmetries: {scaling_dim}")
    print("-" * 30)
    print("Total Dimension =",
          f" {time_translation_dim} + {space_translation_dim} + {rotation_dim} + {galilean_boost_dim} + {scaling_dim} =",
          f" {total_dimension}")

if __name__ == "__main__":
    calculate_lie_group_dimension()