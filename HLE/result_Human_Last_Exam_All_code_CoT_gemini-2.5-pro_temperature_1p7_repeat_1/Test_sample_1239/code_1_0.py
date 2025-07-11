def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the maximal finite-dimensional Lie group of
    symmetries for the 3D incompressible Navier-Stokes equations.

    The dimension is the sum of the number of generators for each type of
    symmetry transformation that leaves the equations invariant.
    """

    # The incompressible Navier-Stokes equations in R^3 are:
    # ∂u/∂t + (u ⋅ ∇)u = -∇p + ν∇²u
    # ∇ ⋅ u = 0

    print("The dimension of the Lie group of symmetries is calculated by summing the number of independent generators for each type of symmetry.")
    print("-" * 20)

    # 1. Time translation (t -> t + a)
    time_translation = 1
    print(f"Time translations: {time_translation} generator")

    # 2. Space translations (x -> x + a) in 3 dimensions
    space_translations = 3
    print(f"Space translations (in x, y, z): {space_translations} generators")

    # 3. Spatial rotations (SO(3) group) in 3 dimensions
    spatial_rotations = 3
    print(f"Spatial rotations (about x, y, z axes): {spatial_rotations} generators")

    # 4. Galilean boosts (transforming to a moving frame) in 3 dimensions
    galilean_boosts = 3
    print(f"Galilean boosts (in x, y, z directions): {galilean_boosts} generators")

    # 5. Scaling transformation (t, x, u, p are scaled)
    scaling_transformation = 1
    print(f"Scaling transformations: {scaling_transformation} generator")
    
    # 6. Pressure shift (p -> p + c), since only ∇p appears in the equations.
    pressure_shift = 1
    print(f"Pressure shifts: {pressure_shift} generator")

    print("-" * 20)

    # Calculate the total dimension
    total_dimension = (time_translation + space_translations + spatial_rotations +
                       galilean_boosts + scaling_transformation + pressure_shift)

    # Print the final equation as requested
    print("The final calculation is the sum of these dimensions:")
    print(f"{time_translation} + {space_translations} + {spatial_rotations} + {galilean_boosts} + {scaling_transformation} + {pressure_shift} = {total_dimension}")
    
    print("\nTherefore, the dimension of the Lie group of symmetries is:")
    print(total_dimension)


calculate_navier_stokes_symmetry_dimension()