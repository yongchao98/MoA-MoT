import sympy

def solve_fluid_interface_equation():
    """
    This function derives and displays the coefficients A(r) and B(r)
    for the governing linear equation of the fluid interface shape ξ(r).
    """
    # Define the symbolic variables used in the equation.
    # r: radial position
    # gamma: surface tension
    # xi(r): fluid interface displacement as a function of r
    r = sympy.Symbol('r', positive=True)
    gamma = sympy.Symbol('gamma')
    xi = sympy.Function('xi')(r)
    
    # Based on the linearization of the Young-Laplace equation in cylindrical
    # coordinates, the governing differential equation is derived as:
    # gamma * (d^2ξ/dr^2) + (gamma/r) * (dξ/dr) - P_elec(r) = 0
    # where P_elec(r) is the pressure induced by the electric field.

    # We match this to the form:
    # A(r) * (d^2ξ/dr^2) + B(r) * (dξ/dr) + C(r, ξ(r)) = 0

    # From this comparison, we identify the coefficients A(r) and B(r).
    A_r = gamma
    B_r = gamma / r
    
    # The problem asks for A(r) and B(r). We will print them.
    # The instruction to "output each number in the final equation" is
    # interpreted as outputting each component/coefficient.

    print("The governing linear equation for the interfacial shape ξ(r) is of the form:")
    print("A(r) * ξ''(r) + B(r) * ξ'(r) + C(r, ξ) = 0\n")
    
    print("--- Coefficient A(r) ---")
    print("This is the coefficient of the second derivative term, ξ''(r).")
    print(f"A(r) = {A_r}\n")

    print("--- Coefficient B(r) ---")
    print("This is the coefficient of the first derivative term, ξ'(r).")
    print(f"B(r) = {B_r}\n")

    # For completeness, let's also show the structure of the C term and the full equation.
    P_elec = sympy.Function('P_elec')(r)
    C_r_xi = -P_elec
    
    print("--- Term C(r, ξ) ---")
    print("This term includes external forces, in this case, the electrostatic pressure.")
    print(f"C(r, ξ) = {C_r_xi}\n")
    
    print("--- Full Governing Equation ---")
    # We construct and print the full equation for clarity.
    equation = sympy.Eq(A_r * xi.diff(r, 2) + B_r * xi.diff(r, 1) + C_r_xi, 0)
    sympy.pprint(equation, use_unicode=True)

if __name__ == "__main__":
    solve_fluid_interface_equation()