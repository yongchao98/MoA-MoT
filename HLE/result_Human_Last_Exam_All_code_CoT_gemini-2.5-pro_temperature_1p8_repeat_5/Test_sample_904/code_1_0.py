import sympy as sp

def solve_fluid_interface_equation():
    """
    This function derives and prints the coefficients A(r) and B(r)
    for the governing differential equation of the fluid interface.
    """
    # Define the symbols used in the problem
    r = sp.Symbol('r')        # Radial position
    gamma = sp.Symbol('γ')    # Surface tension
    xi = sp.Function('ξ')(r)  # Interface displacement

    # Based on the Young-Laplace equation for an axisymmetric interface under the
    # small slope approximation, the pressure balance is:
    # γ * (d²ξ/dr² + (1/r) * dξ/dr) + C(r, ξ) = 0
    # where C(r, ξ) represents the electrostatic pressure term.

    # We compare this to the general form:
    # A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ) = 0

    # From this comparison, we can identify the coefficients A(r) and B(r).
    # A(r) is the coefficient of the second derivative term.
    A_r = gamma

    # B(r) is the coefficient of the first derivative term.
    B_r = gamma / r
    
    # Print the results in a clear format
    print("The derivation is based on the linearized Young-Laplace equation for an axisymmetric fluid interface.")
    print("The governing differential equation is of the form: A(r) * ξ''(r) + B(r) * ξ'(r) + C(r, ξ) = 0")
    print("-" * 50)
    
    print("The coefficient A(r) is the term multiplying the second derivative d²ξ/dr².")
    print(f"A(r) = {A_r}")
    
    print("\nThe coefficient B(r) is the term multiplying the first derivative dξ/dr.")
    print(f"B(r) = {B_r}")

    print("\nTherefore, the final equation takes the form:")
    # Using unicode for better representation
    print(f"{A_r} * d²ξ/dr² + ({B_r}) * dξ/dr + C(r, ξ) = 0")


solve_fluid_interface_equation()