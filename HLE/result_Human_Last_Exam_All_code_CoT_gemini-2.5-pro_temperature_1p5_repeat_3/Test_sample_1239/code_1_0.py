def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the maximal finite-dimensional Lie group
    of symmetries for the 3D incompressible Navier-Stokes equations.
    """

    # Dictionary mapping each type of symmetry to its dimension
    symmetries = {
        "Time translation": 1,
        "Space translations": 3,
        "Spatial rotations (in 3D)": 3,
        "Galilean boosts": 3,
        "Scaling (Dilation)": 1,
        "Constant pressure shift": 1
    }

    print("The dimension of the Lie group of symmetries is determined by summing the number of independent generators for each type of transformation:")
    print("-" * 60)

    total_dimension = 0
    equation_parts = []

    for name, dim in symmetries.items():
        print(f"- {name}: {dim} parameter(s)")
        total_dimension += dim
        equation_parts.append(str(dim))

    print("-" * 60)
    # The prompt requires printing the equation with each number explicitly.
    final_equation_str = " + ".join(equation_parts)
    print(f"The total dimension is the sum: {final_equation_str}")
    print(f"Total Dimension = {total_dimension}")


if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()
