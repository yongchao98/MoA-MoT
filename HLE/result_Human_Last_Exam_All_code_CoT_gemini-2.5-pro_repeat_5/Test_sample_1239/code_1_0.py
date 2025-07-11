def calculate_navier_stokes_symmetry_dimension():
    """
    Calculates the dimension of the Lie group of symmetries for the
    3D incompressible Navier-Stokes equations.

    This is done by summing the dimensions of each independent symmetry
    transformation that leaves the equations invariant.
    """

    # Dictionary mapping each symmetry to its dimension (number of parameters)
    symmetries = {
        "Time translation (t -> t + a)": 1,
        "Spatial translations (x -> x + b)": 3,
        "Spatial rotations (x -> Rx)": 3,
        "Galilean boosts (x -> x + vt, u -> u + v)": 3,
        "Scaling (t,x,u,p -> a^2*t, a*x, a^-1*u, a^-2*p)": 1,
        "Pressure shift (p -> p + c)": 1,
    }

    print("The dimension of the Lie group of symmetries for the 3D incompressible Navier-Stokes equations is calculated by summing the dimensions of its independent symmetry subgroups:")
    print("-" * 80)

    total_dimension = 0
    equation_parts = []

    for name, dim in symmetries.items():
        print(f"- {name}: {dim} dimension(s)")
        total_dimension += dim
        equation_parts.append(str(dim))
    
    print("-" * 80)
    
    # Format the final equation string
    equation_str = " + ".join(equation_parts)
    
    print("The total dimension is the sum of these dimensions:")
    print(f"{equation_str} = {total_dimension}")

# Execute the calculation
calculate_navier_stokes_symmetry_dimension()