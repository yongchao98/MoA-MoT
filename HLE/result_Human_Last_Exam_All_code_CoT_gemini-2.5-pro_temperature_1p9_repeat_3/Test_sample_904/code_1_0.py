import sympy as sp

def solve_fluid_interface_equation():
    """
    This function explains the derivation of the governing linear equation
    for the fluid interface and identifies the coefficients A(r) and B(r).
    """

    # Define symbols
    r, gamma = sp.symbols('r gamma')
    xi = sp.Function('xi')(r)

    print("Derivation of the Governing Equation for the Fluid Interface:")
    print("="*60)

    # Step 1: Start with the linearized capillary pressure equation.
    print("Step 1: The pressure balance equation for the interface is given by the sum of capillary pressure and external (electrostatic) pressure.")
    print("The linearized capillary pressure for an axisymmetric interface is:")
    print(f"ΔP_cap = γ * (d²ξ/dr² + (1/r) * dξ/dr)\n")

    # Step 2: Full governing equation
    print("Step 2: The full pressure balance equation (neglecting gravity) is:")
    print("ΔP_cap + ΔP_ext(r, ξ) = 0")
    print("Substituting the linearized capillary pressure gives:")
    print(f"γ * d²ξ/dr² + (γ/r) * dξ/dr + C(r, ξ) = 0")
    print("where C(r, ξ) represents the external electrostatic pressure term ΔP_ext.\n")

    # Step 3: Identify coefficients
    print("Step 3: Comparing this with the general form A(r)d²ξ/dr² + B(r)dξ/dr + C(r, ξ) = 0, we can identify the coefficients A(r) and B(r).\n")

    # Define A(r) and B(r)
    A_r = gamma
    B_r = gamma / r

    # Print the results
    print("Final Result:")
    print("-" * 20)
    print(f"The coefficient A(r) is the coefficient of the d²ξ/dr² term.")
    print(f"A(r) = {A_r}")
    print("\n")
    print(f"The coefficient B(r) is the coefficient of the dξ/dr term.")
    print(f"B(r) = {B_r}")
    print("-" * 20)

    # Output for the final answer format
    print("\nFinal Answer in requested format:")
    # We are asked to output each number/symbol in the final equation.
    # For A(r) = γ, the symbol is γ.
    # For B(r) = γ/r, the symbols are γ and r.
    final_answer_A = "A(r) = γ"
    final_answer_B = "B(r) = γ / r"
    print(f"The expression for A(r) is: {final_answer_A}")
    print(f"The expression for B(r) is: {final_answer_B}")

if __name__ == '__main__':
    solve_fluid_interface_equation()
