import sys

def calculate_lie_group_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    3D incompressible Navier-Stokes equations.
    """
    # Dimension for time translation (t -> t + a)
    time_translation_dim = 1
    
    # Dimension for spatial translations in R^3 (x -> x + b)
    space_translations_dim = 3
    
    # Dimension for spatial rotations in R^3 (rotations around x, y, z axes)
    space_rotations_dim = 3
    
    # Dimension for Galilean boosts (change to a moving reference frame)
    galilean_boosts_dim = 3
    
    # Dimension for the scaling transformation
    scaling_dim = 1

    # The total dimension is the sum of the dimensions of the independent subgroups
    total_dimension = (
        time_translation_dim +
        space_translations_dim +
        space_rotations_dim +
        galilean_boosts_dim +
        scaling_dim
    )

    print("The dimension is the sum of the dimensions from each symmetry subgroup:")
    print(f"Time Translation: {time_translation_dim}")
    print(f"Spatial Translations: {space_translations_dim}")
    print(f"Spatial Rotations: {space_rotations_dim}")
    print(f"Galilean Boosts: {galilean_boosts_dim}")
    print(f"Scaling Transformation: {scaling_dim}")
    print("\nFinal Calculation:")
    print(
        f"{time_translation_dim} + {space_translations_dim} + {space_rotations_dim} + "
        f"{galilean_boosts_dim} + {scaling_dim} = {total_dimension}"
    )

if __name__ == "__main__":
    calculate_lie_group_dimension()