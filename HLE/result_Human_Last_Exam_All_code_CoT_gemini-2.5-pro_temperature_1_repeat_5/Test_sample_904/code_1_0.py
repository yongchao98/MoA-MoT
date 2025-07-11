import sympy as sp

def solve_fluid_equation():
    """
    This function uses the sympy library to symbolically represent the derivation
    of the governing equation for the fluid interface and prints the coefficients.
    """
    # Define symbols and functions for the derivation
    r = sp.symbols('r', real=True, positive=True)
    
    # gamma is the surface tension, which depends on the radial position r
    # through its dependence on the electric potential, V(r).
    gamma = sp.Function('gamma')(r)
    
    # xi is the displacement of the fluid interface
    xi = sp.Function('xi')(r)
    
    # E is the magnitude of the electric field
    E = sp.Function('E')(r)

    # epsilon_a and epsilon_b are the permittivities of the upper and lower fluids
    epsilon_a, epsilon_b = sp.symbols('epsilon_a epsilon_b')

    # The governing equation for the interface shape xi(r) is derived from the
    # principle of least action (minimizing the total energy of the system).
    # The resulting linear differential equation is of the form:
    # A(r) * d^2(xi)/dr^2 + B(r) * d(xi)/dr + C(r) = 0
    #
    # Our derivation yields the following expressions for the coefficients.

    # Coefficient A(r) is the surface tension, which varies with r due to the potential V(r).
    A_r = gamma

    # Coefficient B(r) includes a term from the geometry (gamma/r) and a term from
    # the spatial variation of the surface tension itself.
    B_r = gamma / r + sp.Derivative(gamma, r)

    # For completeness, the term C(r) represents the electrostatic pressure on the interface.
    C_r = - (1/2) * (epsilon_a - epsilon_b) * E**2

    # The final equation is A_r * xi'' + B_r * xi' + C_r = 0
    
    print("The governing linear equation for the interfacial shape ξ(r) is:")
    print("A(r) * d²ξ/dr² + B(r) * dξ/dr + C(r, ξ(r)) = 0\n")
    print("Based on the derivation, the expressions for the coefficients A(r) and B(r) are:\n")

    # Print the expression for A(r)
    print("A(r) = " + str(A_r))
    print("     (This is the surface tension, which depends on r via the electric potential V(r))\n")

    # Print the expression for B(r)
    print("B(r) = " + str(B_r))
    print("     (This includes the geometric term gamma(r)/r and a term from the spatial change in surface tension)")

    # The problem asks for A(r) and B(r). Here they are.
    # The prompt also asks to "output each number in the final equation!". To fulfill this, we show C(r) as well.
    print("\nFor completeness, the expression for C(r, ξ(r)) is:")
    print("C(r, ξ(r)) = " + str(C_r))
    print("             (This term represents the electrostatic pressure on the interface and is independent of ξ in this linear model)")

if __name__ == '__main__':
    solve_fluid_equation()