def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the incompressible Navier-Stokes equations in R^3.
    """

    print("The full Lie group of symmetries of the Navier-Stokes equations is infinite-dimensional.")
    print("This is due to a gauge symmetry allowing the pressure to be shifted by an arbitrary function of time.")
    print("However, the standard question refers to the dimension of the maximal finite-dimensional symmetry subgroup.")
    print("\nThis finite-dimensional group is composed of the following physical symmetries:")

    # Dimension from time translation: t -> t + a
    dim_time_translation = 1
    print(f"- Time translation: {dim_time_translation} dimension")

    # Dimension from space translations: x -> x + a (vector in R^3)
    dim_space_translation = 3
    print(f"- Space translations: {dim_space_translation} dimensions")

    # Dimension from space rotations: x -> Rx (R in SO(3))
    dim_space_rotation = 3
    print(f"- Space rotations: {dim_space_rotation} dimensions")

    # Dimension from Galilean boosts: x -> x + vt, u -> u + v
    dim_galilean_boosts = 3
    print(f"- Galilean boosts: {dim_galilean_boosts} dimensions")

    # Dimension from scaling transformation (dilation)
    dim_scaling = 1
    print(f"- Scaling transformation: {dim_scaling} dimension")
    
    # Dimension from constant pressure shift: p -> p + c
    dim_pressure_shift = 1
    print(f"- Constant pressure shift: {dim_pressure_shift} dimension")

    # Calculate the total dimension
    total_dimension = (
        dim_time_translation +
        dim_space_translation +
        dim_space_rotation +
        dim_galilean_boosts +
        dim_scaling +
        dim_pressure_shift
    )
    
    # Create the equation string
    equation = (
        f"{dim_time_translation} + {dim_space_translation} + "
        f"{dim_space_rotation} + {dim_galilean_boosts} + "
        f"{dim_scaling} + {dim_pressure_shift}"
    )

    print(f"\nThe total dimension is the sum of these individual dimensions:")
    print(f"Total dimension = {equation} = {total_dimension}")

if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()