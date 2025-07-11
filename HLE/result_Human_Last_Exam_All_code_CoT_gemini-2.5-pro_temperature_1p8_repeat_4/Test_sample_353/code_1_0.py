import sympy

def solve_bdf4_stability_angle():
    """
    Calculates the exact A(alpha)-stability angle for the BDF4 method.
    """
    # Step 1: Define the characteristic polynomials for BDF4.
    # The BDF4 method is given by the formula:
    # 25/12 * y_{n+4} - 4 * y_{n+3} + 3 * y_{n+2} - 4/3 * y_{n+1} + 1/4 * y_n = h * f_{n+4}
    # To work with integer coefficients, we multiply by 12, leading to:
    # 25*y_{n+4} - 48*y_{n+3} + 36*y_{n+2} - 16*y_{n+1} + 3*y_n = 12*h*f_{n+4}
    # From this, we derive the characteristic polynomials rho(s) and sigma(s).
    s = sympy.Symbol('s')
    rho = 25 * s**4 - 48 * s**3 + 36 * s**2 - 16 * s + 3
    sigma = 12 * s**4
    
    print("Step 1: The characteristic polynomials for BDF4 are:")
    print(f"rho(s) = {rho}")
    print(f"sigma(s) = {sigma}\n")

    # Step 2: Formulate the equation for extremal angles.
    # The condition for the argument of z(theta) to be extremal leads to a
    # polynomial equation for x = cos(theta). The derivation yields:
    x = sympy.Symbol('x')
    poly_eqn = 5*x**4 - 16*x**3 + 18*x**2 - 8*x + 1
    
    print("Step 2: The polynomial equation for x = cos(theta) at extremal angles is:")
    print(f"{poly_eqn} = 0\n")

    # Step 3: Solve the polynomial for x = cos(theta).
    roots = sympy.solve(poly_eqn, x)
    # The root x=1 (cos(theta)=1) corresponds to theta=0, where z(0)=0. This is a trivial case.
    # The other root determines the angle alpha.
    cos_theta_val = [r for r in roots if r != 1][0]
    
    print("Step 3: The roots of the polynomial are [1, 1, 1, 1/5].")
    print(f"The non-trivial root is cos(theta) = {cos_theta_val}\n")

    # Step 4: Calculate the complex value z at the extremal angle.
    sin_theta_val = sympy.sqrt(1 - cos_theta_val**2)
    
    # We use z(theta) = rho(exp(i*theta)) / sigma(exp(i*theta)) which can be expressed as:
    # z(theta) = (1/12) * (25 - 48*e^{-it} + 36*e^{-2it} - 16*e^{-3it} + 3*e^{-4it})
    # We calculate the real and imaginary parts using cos(k*theta) and sin(k*theta).
    cos_t = cos_theta_val
    sin_t = sin_theta_val
    cos_2t = 2 * cos_t**2 - 1
    sin_2t = 2 * sin_t * cos_t
    cos_3t = 4 * cos_t**3 - 3 * cos_t
    sin_3t = sin_t * (3 - 4*sin_t**2)
    cos_4t = 2 * cos_2t**2 - 1
    sin_4t = 2 * sin_2t * cos_2t
    
    re_z = (1/sympy.S(12)) * (25 - 48*cos_t + 36*cos_2t - 16*cos_3t + 3*cos_4t)
    im_z = (1/sympy.S(12)) * (48*sin_t - 36*sin_2t + 16*sin_3t - 3*sin_4t)
    
    # Simplify the expressions
    re_z_simplified = sympy.simplify(re_z)
    im_z_simplified = sympy.simplify(im_z)
    
    print("Step 4: The complex number z at this angle is:")
    print(f"Re(z) = {re_z_simplified}")
    print(f"Im(z) = {im_z_simplified}\n")
    
    # Step 5: Calculate the stability angle alpha.
    # The angle alpha is given by alpha = arctan(-Im(z)/Re(z)).
    tan_alpha = -im_z_simplified / re_z_simplified
    tan_alpha_simplified = sympy.cancel(tan_alpha)
    
    print("Step 5: The stability angle alpha is given by arctan(-Im(z)/Re(z)).")
    print(f"tan(alpha) = {tan_alpha_simplified}\n")
    
    # Final step: Express the angle alpha using arctan and output the numbers in the final expression.
    print("--------------------------------------------------")
    print("Final Answer: The exact value of the angle alpha")
    print("--------------------------------------------------")
    
    # Extract the coefficients for clear printing
    num, den = tan_alpha_simplified.as_numer_denom()
    num_coeff = num / sympy.sqrt(6) # The numerator has a sqrt(6) factor
    
    # Print the equation with its numerical components
    print("The final equation for alpha is:")
    print(f"alpha = arctan( (A * sqrt(B)) / C )")
    print("\nWhere the numerical values are:")
    print(f"A = {num_coeff}")
    print(f"B = 6")
    print(f"C = {den}")

if __name__ == "__main__":
    solve_bdf4_stability_angle()