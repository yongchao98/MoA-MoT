def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    incompressible Navier-Stokes equations in R^3.
    
    The total dimension is the sum of the dimensions of the individual
    symmetry transformations.
    """

    # 1. Time translations (t -> t + a)
    time_translations_dim = 1
    print(f"Dimension of time translations: {time_translations_dim}")

    # 2. Space translations (x -> x + b) in 3D
    space_translations_dim = 3
    print(f"Dimension of space translations: {space_translations_dim}")

    # 3. Rotations (SO(3) group) in 3D
    rotations_dim = 3
    print(f"Dimension of spatial rotations: {rotations_dim}")

    # 4. Galilean boosts (changing inertial frame) in 3D
    galilean_boosts_dim = 3
    print(f"Dimension of Galilean boosts: {galilean_boosts_dim}")

    # 5. Scaling transformations
    scaling_dim = 1
    print(f"Dimension of scaling transformations: {scaling_dim}")

    # The total dimension is the sum of these individual dimensions.
    total_dimension = (time_translations_dim + 
                       space_translations_dim + 
                       rotations_dim + 
                       galilean_boosts_dim + 
                       scaling_dim)

    # Print the final equation with each number
    print("\n--- Total Dimension Calculation ---")
    print(f"Total Dimension = {time_translations_dim} (time) + {space_translations_dim} (space) + {rotations_dim} (rotations) + {galilean_boosts_dim} (boosts) + {scaling_dim} (scaling)")
    print(f"Total Dimension = {total_dimension}")

if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()
