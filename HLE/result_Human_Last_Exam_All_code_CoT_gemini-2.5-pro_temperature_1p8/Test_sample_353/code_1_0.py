import numpy as np

def solve_bdf4_stability_angle():
    """
    Calculates the A(alpha)-stability angle for the BDF4 method.
    """
    # The condition for the tangent to the stability boundary passing through the
    # origin can be formulated as a polynomial equation in c = cos(theta).
    #
    # The equation is derived from Re( (xi*rho'(xi)/rho(xi)) - 4 ) = 0 for xi = exp(i*theta).
    # This leads to a high-degree polynomial in c = cos(theta). After simplification
    # and factoring out the trivial solution at c=1 (theta=0), we get the following
    # polynomial of degree 6. We need to find its root c0 between -1 and 1.
    
    # The coefficients of the polynomial equation:
    # 1280*c^6 - 2304*c^5 - 768*c^4 + 1728*c^3 - 124*c^2 - 392*c + 10 = 0
    coeffs = [1280, -2304, -768, 1728, -124, -392, 10]
    roots = np.roots(coeffs)

    # We filter for the real root in the interval (-1, 1).
    c0 = 0
    for r in roots:
        if abs(r.imag) < 1e-9 and -1 < r.real < 1:
            c0 = r.real
            break
            
    # With c0 found, calculate theta0, then the point z0 on the boundary.
    theta0 = np.arccos(c0)
    xi0 = np.exp(1j * theta0)

    # BDF4 rho polynomial coefficients, from xi^4 to xi^0
    rho_coeffs = [25/12, -4, 3, -4/3, 1/4]
    rho_poly = np.poly1d(rho_coeffs)

    z0 = rho_poly(xi0) / (xi0**4)
    
    # The stability angle alpha is given by pi - |arg(z0)|.
    # The problem asks for the answer in terms of arctan.
    # We have tan(alpha) = tan(pi - |arg(z0)|) = -tan(|arg(z0)|).
    # Since z0 is in the third quadrant, arg(z0) is negative, so |arg(z0)|=-arg(z0).
    # tan(alpha) = -tan(-arg(z0)) = tan(arg(z0)) = Im(z0)/Re(z0).
    # Wait, my logic was circular. The simplest formula is tan(alpha) = -Im(z0)/Re(z0)
    # since alpha = pi - angle, where angle is from the positive x-axis.
    tan_alpha = -z0.imag / z0.real
    
    print("The final equation to solve for c = cos(theta) is:")
    print("1280*c^6 - 2304*c^5 - 768*c^4 + 1728*c^3 - 124*c^2 - 392*c + 10 = 0")
    print("\nSolving this numerically, we find the relevant root c0.")
    print(f"The value of c0 is approximately: {c0:.6f}")
    
    print("\nFrom c0, we calculate the tangent of the stability angle alpha.")
    print("The exact value is alpha = arctan(T), where T is:")
    print(f"T = {tan_alpha}")
    print("\nSo the final expression is:")
    print(f"alpha = arctan({tan_alpha:.10f})")

solve_bdf4_stability_angle()
<<<arctan(3.3111322194)>>>