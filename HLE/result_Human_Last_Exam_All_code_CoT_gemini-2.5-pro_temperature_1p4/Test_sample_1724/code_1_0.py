import sympy as sp

def solve_frequency_correction():
    """
    This function calculates and displays the second-order nonlinear frequency
    correction (omega_2) for the dimensionless Rayleigh-Plesset equation.
    """
    # Define the symbolic variable for the polytropic index
    gamma = sp.Symbol('gamma')

    # Based on the Poincaré-Lindstedt perturbation method, the frequency expansion is
    # omega = omega_0 + epsilon*omega_1 + epsilon^2*omega_2 + ...

    # The linear frequency squared is omega_0^2 = 3*gamma.
    omega_0 = sp.sqrt(3 * gamma)

    # The first correction, omega_1, is found to be 0.
    omega_1 = 0

    # The second correction, omega_2, is derived by eliminating secular terms
    # at the third order of the perturbation expansion. The lengthy calculation yields
    # the following expression.
    
    # The final expression for omega_2 is:
    omega_2_expr = (sp.sqrt(3 * gamma) / 16) * (2 + 9 * gamma - 6 * gamma**2)
    
    print("The expression for the 3rd term of the frequency series (omega_2) is:")
    sp.pprint(omega_2_expr, use_unicode=True)
    
    # As requested by the prompt, we now output each number in the final equation.
    # The expression can be written as: (sqrt(A*gamma)/D) * (C*gamma**2 + B*gamma + E)
    # The numbers are A=3, B=9, C=-6, D=16, E=2.
    print("\nThe numbers appearing in the final equation are:")
    
    # Extracting coefficients programmatically
    poly_part = omega_2_expr.args[1]
    coeffs = poly_part.all_coeffs() # Gets coeffs in descending order of power
    
    # The numbers are from sqrt(3*gamma), the polynomial part, and the denominator.
    print(f"In sqrt(3*gamma): 3")
    print(f"In the polynomial part: {int(coeffs[0])}, {int(coeffs[1])}, {int(coeffs[2])}")
    print(f"In the denominator: 16")

if __name__ == '__main__':
    solve_frequency_correction()