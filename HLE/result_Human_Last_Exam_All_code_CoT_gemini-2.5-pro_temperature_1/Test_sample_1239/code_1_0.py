def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the incompressible Navier-Stokes equations in R^3.
    """

    # The Lie group of symmetries is composed of several subgroups. We sum their dimensions.

    # 1. Time translation: t' = t + a
    # This transformation is described by a single parameter 'a'.
    time_translation_dim = 1
    
    # 2. Spatial translations: x' = x + b
    # In R^3, the translation vector 'b' has 3 components.
    spatial_translations_dim = 3
    
    # 3. Spatial rotations: x' = R * x
    # The group of rotations in 3D, SO(3), has 3 parameters (e.g., Euler angles).
    spatial_rotations_dim = 3
    
    # 4. Galilean boosts: x' = x + v*t, u' = u + v
    # The boost velocity vector 'v' has 3 components in R^3.
    galilean_boosts_dim = 3
    
    # 5. Scaling transformation: t' = λ^2*t, x' = λ*x, u' = λ^-1*u
    # This is a one-parameter group defined by the scaling factor 'λ'.
    scaling_dim = 1

    # 6. Constant pressure shift: p' = p + c
    # The pressure only appears as a gradient, so adding a constant 'c' has no effect.
    # This adds one more dimension.
    pressure_shift_dim = 1

    # The total dimension is the sum of the dimensions of these independent transformations.
    total_dimension = (time_translation_dim +
                       spatial_translations_dim +
                       spatial_rotations_dim +
                       galilean_boosts_dim +
                       scaling_dim +
                       pressure_shift_dim)

    print("The total dimension of the Lie group is the sum of the dimensions of its constituent symmetries:")
    print(f"Dimension from time translation      = {time_translation_dim}")
    print(f"Dimension from spatial translations = {spatial_translations_dim}")
    print(f"Dimension from spatial rotations    = {spatial_rotations_dim}")
    print(f"Dimension from Galilean boosts      = {galilean_boosts_dim}")
    print(f"Dimension from scaling              = {scaling_dim}")
    print(f"Dimension from pressure shift       = {pressure_shift_dim}")
    print("-" * 40)
    
    print("Final Calculation:")
    print(f"{time_translation_dim} + {spatial_translations_dim} + {spatial_rotations_dim} + {galilean_boosts_dim} + {scaling_dim} + {pressure_shift_dim} = {total_dimension}")
    print(f"\nThe dimension of the Lie group of symmetries is {total_dimension}.")


if __name__ == '__main__':
    calculate_navier_stokes_symmetry_dimension()