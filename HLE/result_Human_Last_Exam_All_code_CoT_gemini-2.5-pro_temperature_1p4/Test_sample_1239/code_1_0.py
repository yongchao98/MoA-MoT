def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    for the 3D incompressible Navier-Stokes equations.
    """
    # The symmetry group is known as the Galilean-similitude group.
    # We define a dictionary mapping each symmetry type to its dimension.
    symmetries = {
        "Time translation": 1,
        "Spatial translations": 3,
        "Spatial rotations": 3,
        "Galilean boosts": 3,
        "Scaling (Dilation)": 1
    }

    print("The finite-dimensional Lie group of symmetries for the incompressible Navier-Stokes equations in R^3 is composed of the following transformations:")
    
    dimensions = []
    for name, dim in symmetries.items():
        print(f"- {name}: Contributes {dim} dimension(s).")
        dimensions.append(dim)
        
    total_dimension = sum(dimensions)
    
    # Constructing the equation string as requested.
    # The map(str, dimensions) converts each number in the list to a string.
    # " + ".join(...) joins these strings with " + " in between.
    equation_string = " + ".join(map(str, dimensions))
    
    print("\nThe total dimension is the sum of the dimensions of these individual symmetries.")
    print("The calculation is:")
    print(f"{equation_string} = {total_dimension}")
    print("\nNote: This is the dimension of the maximal finite-dimensional part of the symmetry group. An additional infinite-dimensional symmetry related to pressure shifts also exists.")

if __name__ == "__main__":
    calculate_navier_stokes_symmetry_dimension()