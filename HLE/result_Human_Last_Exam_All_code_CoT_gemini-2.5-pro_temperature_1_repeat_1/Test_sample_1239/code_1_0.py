def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries of the
    3D incompressible Navier-Stokes equations.
    """

    # The dimension of a Lie group of symmetries is the number of independent
    # parameters needed to describe the transformations. We sum the dimensions
    # corresponding to each type of symmetry.

    # 1. Time translation: t -> t + a
    # This transformation is described by one parameter.
    time_translation_dim = 1

    # 2. Space translations: x -> x + b
    # In R^3, this transformation is described by a 3D vector, so 3 parameters.
    space_translations_dim = 3

    # 3. Space rotations: x -> R*x (where R is a rotation matrix)
    # The group of rotations in 3D (SO(3)) has 3 dimensions (e.g., Euler angles).
    space_rotations_dim = 3

    # 4. Galilean boosts: x -> x + v*t, u -> u + v
    # This transformation to a moving frame is described by a 3D velocity vector, so 3 parameters.
    galilean_boosts_dim = 3

    # 5. Scaling (Dilation): t -> a^2*t, x -> a*x, u -> a^-1*u
    # This transformation is described by one scaling parameter 'a'.
    scaling_dim = 1

    # The total dimension is the sum of the dimensions of these individual symmetries.
    total_dimension = (time_translation_dim +
                       space_translations_dim +
                       space_rotations_dim +
                       galilean_boosts_dim +
                       scaling_dim)

    print("The dimension of the Lie group of symmetries for the 3D incompressible Navier-Stokes equations is calculated by summing the dimensions of its component symmetries:")
    print(f"   Time Translation: {time_translation_dim}")
    print(f"   Space Translations: {space_translations_dim}")
    print(f"   Space Rotations: {space_rotations_dim}")
    print(f"   Galilean Boosts: {galilean_boosts_dim}")
    print(f"   Scaling: {scaling_dim}")
    print("-" * 40)
    print("Final Equation:")
    print(f"   Total Dimension = {time_translation_dim} + {space_translations_dim} + {space_rotations_dim} + {galilean_boosts_dim} + {scaling_dim}")
    print(f"   Total Dimension = {total_dimension}")

calculate_navier_stokes_symmetry_dimension()