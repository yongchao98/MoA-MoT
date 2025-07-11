import sys

def solve_lie_group_dimension():
    """
    Calculates and explains the dimension of the finite-dimensional Lie group
    of symmetries of the 3D incompressible Navier-Stokes equations.
    """

    # The finite-dimensional symmetries and their respective dimensions (number of generators).
    # This is a classical result from the symmetry analysis of differential equations.
    symmetries = {
        "Time translation": 1,
        "Space translations": 3,
        "Spatial rotations": 3,
        "Galilean boosts": 3,
        "Scaling (dilations)": 1
    }

    print("The dimension of the maximal finite-dimensional Lie group of symmetries for the")
    print("3D incompressible Navier-Stokes equations is the sum of the dimensions")
    print("of its component symmetry groups:")
    print("-" * 50)
    for name, dim in symmetries.items():
        print(f"Dimension from {name}: {dim}")
    print("-" * 50)

    # Build and print the final equation
    dimensions = list(symmetries.values())
    total_dimension = sum(dimensions)
    
    equation_parts = [str(d) for d in dimensions]
    equation_string = " + ".join(equation_parts)
    
    print("The total dimension is calculated by the sum:")
    print(f"{equation_string} = {total_dimension}")

if __name__ == "__main__":
    solve_lie_group_dimension()
    # The final answer is 11.
    # We add this special format for the final answer as requested.
    # To avoid printing the marker to stdout, we write to stderr.
    sys.stderr.write("<<<11>>>")