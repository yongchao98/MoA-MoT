def solve_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the 3D incompressible Navier-Stokes equations.
    """
    symmetries = {
        "Time translation": 1,
        "Space translations": 3,
        "Spatial rotations": 3,
        "Galilean boosts": 3,
        "Scaling (Dilation)": 1,
        "Constant pressure shift": 1
    }

    print("The dimension of the Lie group of symmetries of the 3D incompressible Navier-Stokes equations is determined by summing the dimensions of its independent continuous symmetry transformations.")
    print("\nThese transformations and their respective dimensions are:")

    total_dimension = 0
    dimension_parts = []
    
    for name, dim in symmetries.items():
        print(f"- {name}: {dim} dimension(s)")
        total_dimension += dim
        dimension_parts.append(str(dim))

    equation_str = " + ".join(dimension_parts)
    print("\nThe total dimension is the sum of these individual dimensions:")
    print(f"Total Dimension = {equation_str} = {total_dimension}")

if __name__ == "__main__":
    solve_navier_stokes_symmetry_dimension()
