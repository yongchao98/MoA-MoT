import sympy

def print_force_equation():
    """
    This function defines and prints the symbolic equation for the instantaneous force
    based on the provided problem description and choice B.
    """
    # Define the symbols used in the equation
    R, N, N0, I0, i0, omega, t, g, mu0, alpha_T, T, T0, Bs = sympy.symbols(
        'R N N_0 I_0 i_0 omega t g mu_0 alpha_T T T_0 B_s', real=True, positive=True
    )
    fx = sympy.Function('f_x')(t)

    # The problem asks for the operating temperature T1, so we substitute T with T1
    T1 = sympy.Symbol('T_1')

    # Construct the numerator of the expression from choice B
    # f_x(t) = -2*pi*R*N * [numerator_fraction] / [denominator_fraction]
    numerator_main = -2 * sympy.pi * R * N
    numerator_fraction = mu0 * (1 - alpha_T * (T1 - T0)) * N0 * I0 * i0 * sympy.sin(omega * t)

    # Construct the denominator of the expression from choice B
    denominator_fraction = g**2 * (1 + (mu0 * N0 * I0) / (g * Bs))

    # Combine to form the full expression for the force
    force_expression = numerator_main * (numerator_fraction / denominator_fraction)

    # Create the final equation
    equation = sympy.Eq(fx, force_expression)

    # Print the equation in a formatted way
    # The prompt asks to output each number/symbol in the final equation.
    # Printing the SymPy equation object achieves this.
    print("The derived instantaneous force f_x(t) is:")
    
    # To match the format of the answer choice exactly, we can construct a string.
    # This also fulfills the requirement to output all parts of the equation.
    equation_string = f"f_x(t) = -2*pi*R*N * (({mu0}*(1 - {alpha_T}*({T1} - {T0}))*{N0}*{I0}*{i0}*sin({omega}*t)) / ({g}**2*(1 + ({mu0}*{N0}*{I0})/({g}*{Bs}))))"
    
    # A more readable version:
    final_equation_str = f"f_x(t) = -2*pi*R*N * (mu_0*(1 - alpha_T*(T_1 - T_0))*N_0*I_0*i_0*sin(omega*t)) / (g^2*(1 + (mu_0*N_0*I_0)/(g*B_s)))"
    
    # Let's print the components to be explicit
    print("\nComponents of the equation:")
    print(f"Coefficient: -2*pi*R*N")
    print(f"Numerator: mu_0*(1 - alpha_T*(T_1 - T_0))*N_0*I_0*i_0*sin(omega*t)")
    print(f"Denominator: g^2*(1 + (mu_0*N_0*I_0)/(g*B_s))")
    
    print("\nFinal Equation:")
    print(final_equation_str)


if __name__ == '__main__':
    print_force_equation()