import collections

def solve_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries of the 3D
    incompressible Navier-Stokes equations.

    The dimension is the sum of the number of generators for each
    independent symmetry transformation.
    """

    # Using an OrderedDict to maintain the order for printing
    symmetries = collections.OrderedDict([
        ("Time translation", 1),
        ("Spatial translations", 3),
        ("Spatial rotations", 3),
        ("Galilean boosts", 3),
        ("Scaling (dilation)", 1),
        ("Constant pressure shift", 1)
    ])

    print("The dimension of the Lie group of symmetries for the 3D incompressible Navier-Stokes equations is calculated by summing the dimensions of its constituent symmetry generators:")
    print("-" * 80)

    for name, dim in symmetries.items():
        print(f"{name:<25}: {dim} dimension(s)")

    print("-" * 80)

    # Calculate the total dimension
    total_dimension = sum(symmetries.values())

    # Create the equation string
    equation_parts = [str(dim) for dim in symmetries.values()]
    equation_str = " + ".join(equation_parts)

    print("The final calculation is:")
    print(f"Total Dimension = {equation_str} = {total_dimension}")


if __name__ == "__main__":
    solve_navier_stokes_symmetry_dimension()