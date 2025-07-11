import math

def generate_potential_expression():
    """
    This function generates and prints the expression for the EDL potential distribution.
    """

    # Define symbols as strings for readable output
    psi_y = "psi(y)"
    z_1 = "z_1"
    beta = "beta"
    k = "k"
    H = "H"
    y = "y"

    print("Derivation of the Electrical Double-Layer (EDL) Potential Distribution:")
    print("-" * 70)
    print("1. Governing Equation (Linearized Poisson-Boltzmann): d^2(psi)/dy^2 = k^2 * psi")
    print(f"2. Boundary Conditions:")
    print(f"   - At top wall (y = H/2): psi(H/2) = 0 (since z_2=0)")
    print(f"   - At bottom wall (y = -H/2): psi(-H/2) = z_a1 = {z_1}*(1 + {beta}*{k})")
    print("-" * 70)

    # Interpreting "output each number in the final equation" as outputting
    # the functional components of the equation.
    print("The final expression for the electrical potential psi(y) is composed of:")

    # Component 1: Slip-dependent zeta potential at the bottom wall
    component1_str = f"{z_1}*(1 + {beta}*{k})"
    print(f"Component 1 (Effective potential at the bottom wall): {component1_str}")

    # Component 2: A term describing the spatial distribution
    component2_str = f"sinh({k}*({H}/2 - {y}))"
    print(f"Component 2 (Spatial distribution term): {component2_str}")

    # Component 3: A normalization factor
    component3_str = f"sinh({k}*{H})"
    print(f"Component 3 (Normalization denominator): {component3_str}")

    print("-" * 70)

    # Construct and print the final, complete expression
    final_expression = f"{psi_y} = ({component1_str} / {component3_str}) * {component2_str}"

    print("The complete expression for the EDL potential distribution is:")
    print(final_expression)


# Run the function to display the results
generate_potential_expression()

# Final answer in the required format
final_answer = "psi(y) = (z_1*(1 + beta*k)/sinh(k*H)) * sinh(k*(H/2 - y))"
# The line below is for the final answer extraction and is not printed to the user.
# <<<psi(y) = (z_1*(1 + beta*k)/sinh(k*H)) * sinh(k*(H/2 - y))>>>