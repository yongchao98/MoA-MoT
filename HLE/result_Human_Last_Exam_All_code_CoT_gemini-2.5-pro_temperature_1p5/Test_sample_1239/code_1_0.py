def solve_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    for the 3D incompressible Navier-Stokes equations.
    """

    print("The full Lie group of symmetries for the incompressible Navier-Stokes equations is infinite-dimensional.")
    print("This arises from the fact that the pressure `p` can be replaced by `p + f(t)` for any arbitrary")
    print("function `f(t)` of time, as pressure only appears under a spatial gradient in the equations.")
    print("\nHowever, the standard question about the 'dimension' refers to the maximal finite-dimensional")
    print("subalgebra of kinematical symmetries. We will calculate the dimension of this group.\n")

    # The components of the finite-dimensional symmetry algebra and their dimensions
    symmetry_components = {
        "Time translations": 1,
        "Spatial translations": 3,
        "Spatial rotations (SO(3))": 3,
        "Galilean boosts": 3,
        "Scaling transformations (dilations)": 1
    }

    print("The finite-dimensional symmetry algebra is spanned by generators for:")
    for component, dim in symmetry_components.items():
        print(f"- {component}: {dim} generator(s)")

    # Calculate the total dimension
    total_dimension = sum(symmetry_components.values())
    
    # Create the equation string as requested
    equation_parts = [str(dim) for dim in symmetry_components.values()]
    equation_string = " + ".join(equation_parts)

    print("\nThe total dimension is the sum of the number of generators for each component:")
    print(f"Total Dimension = {equation_string} = {total_dimension}")

if __name__ == "__main__":
    solve_navier_stokes_symmetry_dimension()