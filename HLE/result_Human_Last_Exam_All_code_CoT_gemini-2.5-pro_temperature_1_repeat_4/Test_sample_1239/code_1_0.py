def solve_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the 3D incompressible Navier-Stokes equations.
    """
    # A dictionary mapping each symmetry transformation to its dimension (number of parameters)
    symmetries = {
        "Time translation": 1,
        "Space translations": 3,
        "Spatial rotations": 3,
        "Galilean boosts": 3,
        "Scaling (Dilation)": 1,
        "Pressure shift": 1
    }

    print("The total dimension of the Lie group is the sum of the dimensions of its constituent symmetry subgroups.")
    print("The symmetries and their respective dimensions are:")
    for name, dim in symmetries.items():
        print(f"- {name}: {dim}")

    # Calculate the total dimension and create the equation string
    dimensions = list(symmetries.values())
    total_dimension = sum(dimensions)
    
    # Format the numbers into a sum equation as requested
    equation_str = " + ".join(map(str, dimensions))
    final_equation = f"{equation_str} = {total_dimension}"

    print("\nThe dimension of the Lie group is the sum of these individual dimensions:")
    print(final_equation)

solve_navier_stokes_symmetry_dimension()