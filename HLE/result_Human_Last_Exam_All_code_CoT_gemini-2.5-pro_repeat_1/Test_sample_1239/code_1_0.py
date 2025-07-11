def solve_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the standard finite-dimensional Lie group of symmetries
    for the 3D incompressible Navier-Stokes equations.

    The full point symmetry group is infinite-dimensional due to the pressure symmetry
    p -> p + g(t). However, the physically non-trivial subgroup, which acts on
    the velocity field, spacetime, and includes a constant pressure shift, is finite-dimensional.
    This script calculates the dimension of this standard subgroup.
    """

    # Define the number of generators for each type of symmetry.
    symmetries = {
        "Time translation": 1,
        "Space translations": 3,
        "Space rotations": 3,
        "Galilean boosts": 3,
        "Scaling transformation": 1,
        "Constant pressure shift": 1
    }

    # Calculate the total dimension by summing the number of generators.
    total_dimension = sum(symmetries.values())

    # Print the breakdown of symmetries and their contributions.
    print("The dimension of the standard Lie group of symmetries for the incompressible Navier-Stokes equations is the sum of the dimensions of its component symmetry groups:")
    for name, dim in symmetries.items():
        print(f"- {name}: {dim} generator(s)")

    # Construct the final calculation string.
    calculation_string = " + ".join(str(d) for d in symmetries.values())

    # Print the final calculation and the result.
    print("\nThe total dimension is calculated as follows:")
    print(f"Dimension = {calculation_string} = {total_dimension}")

if __name__ == "__main__":
    solve_navier_stokes_symmetry_dimension()