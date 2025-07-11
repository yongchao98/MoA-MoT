import sympy
from sympy import symbols, Function, Eq, Derivative

def solve_fluid_equation():
    """
    Symbolically derives the governing linear equation for the fluid interface
    and identifies the coefficients A(r) and B(r).
    """
    # Define the symbolic variables and functions
    r = symbols('r', positive=True, real=True)         # Radial coordinate
    gamma = symbols('gamma', positive=True, real=True) # Surface tension
    xi = Function('xi')(r)                             # Fluid interface displacement
    
    # The external pressure difference (from the electric field)
    # is a function of r. We represent it as C(r) for generality,
    # noting that in the final equation C(r, xi) = -Delta_P_ext(r)
    Delta_P_ext = Function('Delta_P_ext')(r)

    # The general form of the equation is:
    # A(r) * xi''(r) + B(r) * xi'(r) + C(r, xi) = 0

    # Step 1: Start from the linearized Young-Laplace equation in cylindrical coordinates.
    # The pressure from surface tension balances the external pressure.
    # gamma * (Linearized Curvature) = Delta_P_ext
    # Linearized Curvature = xi''(r) + (1/r) * xi'(r)
    
    lhs = gamma * (Derivative(xi, r, 2) + (1/r) * Derivative(xi, r))
    rhs = Delta_P_ext
    
    # Step 2: Form the full differential equation by moving all terms to one side.
    governing_eq_lhs = lhs - rhs
    governing_eq = Eq(governing_eq_lhs, 0)
    
    # Step 3: Identify coefficients A(r) and B(r) by extracting the
    # coefficients of the second and first derivatives of xi(r).
    
    A_r = governing_eq.lhs.coeff(Derivative(xi, r, 2))
    B_r = governing_eq.lhs.coeff(Derivative(xi, r))

    # Print the derived coefficients which form the final equation
    # A(r) * xi'' + B(r) * xi' - Delta_P_ext(r) = 0
    print("Based on the linearization of the Young-Laplace equation in cylindrical coordinates, the coefficients of the governing differential equation are:")
    print(f"A(r) = {A_r}")
    print(f"B(r) = {B_r}")

solve_fluid_equation()