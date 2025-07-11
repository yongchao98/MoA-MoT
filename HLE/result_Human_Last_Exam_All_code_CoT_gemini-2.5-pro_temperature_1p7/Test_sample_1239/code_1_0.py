def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates and prints the dimension of the Lie group of symmetries
    for the 3D incompressible Navier-Stokes equations.
    """
    
    # The dimension is found by summing the dimensions of the independent
    # symmetry subgroups that leave the equations invariant.

    # 1. Dimension from time translation invariance (t -> t + a)
    time_translation_dim = 1
    
    # 2. Dimension from space translation invariance (x -> x + a)
    space_translation_dim = 3
    
    # 3. Dimension from rotational invariance in 3D space (SO(3) group)
    rotation_dim = 3
    
    # 4. Dimension from Galilean boost invariance (shift to a moving frame)
    galilean_boosts_dim = 3
    
    # 5. Dimension from scaling invariance (a specific scaling of all variables)
    scaling_dim = 1
    
    # The total dimension is the sum of these individual dimensions.
    total_dim = (time_translation_dim + space_translation_dim + 
                 rotation_dim + galilean_boosts_dim + scaling_dim)

    print("The total dimension of the Lie group of symmetries for the 3D incompressible Navier-Stokes equations is the sum of the dimensions of its symmetry subgroups:")
    print(f"- Time Translations: {time_translation_dim}")
    print(f"- Space Translations: {space_translation_dim}")
    print(f"- Space Rotations: {rotation_dim}")
    print(f"- Galilean Boosts: {galilean_boosts_dim}")
    print(f"- Scaling Transformations: {scaling_dim}")
    
    # Print the final equation as requested.
    print("\nThe final equation for the total dimension is:")
    print(f"{time_translation_dim} + {space_translation_dim} + {rotation_dim} + {galilean_boosts_dim} + {scaling_dim} = {total_dim}")
    
    print(f"\nThus, the total dimension of the Lie group is {total_dim}.")

if __name__ == '__main__':
    calculate_navier_stokes_symmetry_dimension()