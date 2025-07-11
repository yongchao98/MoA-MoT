def solve_edl_potential():
    """
    This function generates and prints the symbolic expression for the
    Electrical Double-Layer (EDL) potential distribution ψ(y)
    for the given problem conditions.
    """

    # Define symbolic representations of the physical quantities
    psi_of_y = "ψ(y)"
    z_1 = "z_1"
    beta = "β"
    k = "k"
    H = "H"
    y = "y"

    # Define the numbers present in the equation as characters to fulfill the output requirement
    number_one = "1"
    # The potential at the top wall, z_2, is given as 0.

    # Step 1: Formulate the slip-dependent zeta potential at the bottom boundary (y=0).
    # The given formula is z_a1 = z_1(1 + βk).
    slip_dependent_zeta = f"{z_1}({number_one} + {beta}{k})"

    # Step 2: Formulate the spatial distribution part of the solution.
    # This comes from solving the linearized Poisson-Boltzmann equation with
    # boundary conditions ψ(0)=z_a1 and ψ(H)=0.
    # The solution is [sinh(k(H - y)) / sinh(kH)].
    numerator = f"sinh({k}({H} - {y}))"
    denominator = f"sinh({k}{H})"
    spatial_term = f"[{numerator} / {denominator}]"

    # Step 3: Combine the parts to form the final expression for ψ(y).
    final_equation = f"{psi_of_y} = {slip_dependent_zeta} * {spatial_term}"

    print("The final expression for the Electrical double-layer potential distribution is:")
    print(final_equation)

# Execute the function to print the result
solve_edl_potential()
