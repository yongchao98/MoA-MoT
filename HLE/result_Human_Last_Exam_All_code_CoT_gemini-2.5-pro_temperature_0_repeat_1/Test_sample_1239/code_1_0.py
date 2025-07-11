def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the 3D incompressible Navier-Stokes equations.
    """
    # Dimensions of the various independent symmetry transformations
    time_translation_dim = 1
    space_translation_dim = 3
    space_rotation_dim = 3
    galilean_boost_dim = 3
    scaling_dim = 1

    # The total dimension is the sum of the dimensions of these transformations.
    total_dimension = (time_translation_dim +
                       space_translation_dim +
                       space_rotation_dim +
                       galilean_boost_dim +
                       scaling_dim)

    print("The total dimension of the Lie group of symmetries is the sum of the dimensions of its components:")
    print(f"- Time translations: {time_translation_dim}")
    print(f"- Space translations: {space_translation_dim}")
    print(f"- Space rotations: {space_rotation_dim}")
    print(f"- Galilean boosts: {galilean_boost_dim}")
    print(f"- Scaling transformations: {scaling_dim}")
    print("\nThe final dimension is calculated by summing these individual dimensions:")
    
    # Print the final equation with each number, as requested.
    print(f"{time_translation_dim} + {space_translation_dim} + {space_rotation_dim} + {galilean_boost_dim} + {scaling_dim} = {total_dimension}")

calculate_navier_stokes_symmetry_dimension()