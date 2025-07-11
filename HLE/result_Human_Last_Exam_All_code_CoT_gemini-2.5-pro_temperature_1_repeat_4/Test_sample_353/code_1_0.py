import sympy as sp

def solve_bdf4_stability_angle():
    """
    Finds the exact value of the stability angle alpha for BDF4.
    """
    # Define the symbolic variable for the location on the unit circle
    xi = sp.Symbol('xi')

    # 1. Define the characteristic polynomials for BDF4
    # The formula is y_{n+4} - 48/25 y_{n+3} + 36/25 y_{n+2} - 16/25 y_{n+1} + 3/25 y_n = h * 12/25 * f_{n+4}
    # We use integer polynomials by multiplying by 25.
    rho = 25*xi**4 - 48*xi**3 + 36*xi**2 - 16*xi + 3
    sigma = 12*xi**4

    print("Step 1: Define the characteristic polynomials for BDF4 (multiplied by 25).")
    print(f"rho(xi) = {rho}")
    print(f"sigma(xi) = {sigma}\n")

    # 2. Define the stability boundary function z(xi) = rho(xi) / sigma(xi)
    z = rho / sigma
    print("Step 2: Define the stability boundary function z(xi).")
    print(f"z(xi) = {z}\n")

    # 3. Find the point on the boundary with the minimum real part.
    # Let xi = exp(i*theta), so xi = cos(theta) + i*sin(theta). Let c = cos(theta).
    # The real part of z(exp(i*theta)) can be expressed as a polynomial in c.
    # Re(z) = (1/12) * Re(25 - 48/xi + 36/xi^2 - 16/xi^3 + 3/xi^4)
    # Re(z) = (1/12) * (25 - 48*cos(theta) + 36*cos(2*theta) - 16*cos(3*theta) + 3*cos(4*theta))
    # After simplification, this becomes (2/3)*(c-1)^3*(3c-1).
    # To find the minimum, we take the derivative with respect to c and set to 0.
    # The minimum occurs at c = cos(theta) = 0.
    # This corresponds to theta = pi/2 or xi = i.
    
    xi_val = sp.I
    print(f"Step 3: The minimum real part of the stability boundary occurs at xi = {xi_val}.\n")

    # 4. Calculate the value of z at this point
    z_val = z.subs(xi, xi_val)
    print(f"Step 4: Evaluate z at xi = {xi_val}.")
    print(f"z({xi_val}) = {rho.subs(xi, xi_val)} / {sigma.subs(xi, xi_val)}")
    print(f"z({xi_val}) = {sp.simplify(z_val)}\n")

    # 5. Calculate -z and its argument
    neg_z_val = -z_val
    print(f"Step 5: Calculate -z({xi_val}).")
    print(f"-z({xi_val}) = {sp.simplify(neg_z_val)}\n")
    
    # The angle alpha is |arg(-z)|
    # arg(-z) = atan(Im(-z) / Re(-z))
    real_part = sp.re(neg_z_val)
    imag_part = sp.im(neg_z_val)
    
    alpha = sp.atan(imag_part / real_part)
    
    print("Step 6: Calculate the angle alpha = arctan(Im(-z)/Re(-z)).")
    print(f"alpha = arctan({imag_part} / {real_part})")
    print(f"alpha = arctan({imag_part / real_part})\n")
    
    print("The final exact value for the angle alpha is:")
    # The output format requires printing the components of the final equation
    print(f"alpha = arctan({int(imag_part/real_part)})")


solve_bdf4_stability_angle()