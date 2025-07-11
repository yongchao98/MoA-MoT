import sympy as sp

def solve_fluid_equation():
    """
    Derives and prints the coefficients A(r) and B(r) for the governing equation of the fluid interface.
    """
    
    # Define symbolic variables for radius 'r' and surface tension 'gamma'
    r = sp.Symbol('r', positive=True)  # Radial position
    gamma = sp.Symbol('gamma')         # Surface tension

    # --- Derivation Steps ---
    # The governing equation for the interface shape ξ(r) is derived from the pressure balance
    # across the interface. The pressure created by surface tension (Laplace pressure) must
    # balance the pressure exerted by the electric field.
    #
    # P_surface_tension + P_electrostatic = 0
    #
    # The Laplace pressure (P_surface_tension) for an axisymmetric interface z = ξ(r)
    # is given by P_laplace = γ * (mean curvature). In cylindrical coordinates, the linearized
    # mean curvature (for small slopes dξ/dr) is:
    #
    # Mean Curvature ≈ d²ξ/dr² + (1/r) * dξ/dr
    #
    # So, the Laplace pressure is:
    # P_surface_tension ≈ γ * (d²ξ/dr² + (1/r) * dξ/dr)
    #
    # Substituting this into the pressure balance equation gives:
    # γ * (d²ξ/dr² + (1/r) * dξ/dr) + P_electrostatic = 0
    #
    # Distributing γ gives:
    # γ * d²ξ/dr² + (γ/r) * dξ/dr + P_electrostatic = 0
    #
    # This equation is in the required form:
    # A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ(r)) = 0
    #
    # By comparing the terms, we can identify the coefficients A(r) and B(r).
    # The term C(r, ξ(r)) corresponds to the electrostatic pressure.

    # Identify A(r), the coefficient of the second derivative term d²ξ/dr²
    A_r = gamma

    # Identify B(r), the coefficient of the first derivative term dξ/dr
    B_r = gamma / r

    # Print the derived coefficients
    print("The governing linear equation for the interfacial shape ξ(r) is:")
    print("A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ(r)) = 0")
    print("\nBased on the linearization of the Young-Laplace equation, the coefficients are:")
    print("-" * 30)
    # The problem asks to output each component of the final equation
    print(f"The coefficient A(r) is the coefficient of d²ξ/dr².")
    print(f"A(r) = {A_r}")
    print("\n")
    print(f"The coefficient B(r) is the coefficient of dξ/dr.")
    print(f"B(r) = {B_r}")
    print("-" * 30)


if __name__ == "__main__":
    solve_fluid_equation()