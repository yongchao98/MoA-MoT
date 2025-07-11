def solve_navier_stokes_symmetry_dimension():
    """
    Calculates and explains the dimension of the Lie group of symmetries
    of the incompressible Navier-Stokes equations in R^3.
    """

    # The dimension of a Lie group is the number of independent parameters that define its transformations.
    # We find this by summing the dimensions of the independent symmetry transformations
    # that leave the incompressible Navier-Stokes equations invariant.

    # Dictionary to store the dimension of each symmetry component.
    symmetries = {
        'Time translation': 1,        # t -> t + a
        'Space translations': 3,      # x -> x + b (in R^3)
        'Spatial rotations': 3,       # x -> Rx, R in SO(3)
        'Galilean boosts': 3,         # x -> x + vt, u -> u + v (in R^3)
        'Scaling': 1,                 # A specific scaling of t, x, u, and p
        'Pressure shift': 1           # p -> p + c (as only grad(p) appears)
    }

    print("The dimension of the Lie group of symmetries of the incompressible Navier-Stokes equations in R^3 can be calculated by summing the dimensions of its independent symmetry components:")
    
    total_dimension = 0
    equation_parts = []
    
    for name, dim in symmetries.items():
        print(f"- {name}: {dim} dimension(s)")
        total_dimension += dim
        equation_parts.append(str(dim))

    equation_str = " + ".join(equation_parts)
    
    print("\nThe total dimension is the sum of these individual dimensions:")
    print(f"Total dimension = {equation_str} = {total_dimension}")

if __name__ == '__main__':
    solve_navier_stokes_symmetry_dimension()