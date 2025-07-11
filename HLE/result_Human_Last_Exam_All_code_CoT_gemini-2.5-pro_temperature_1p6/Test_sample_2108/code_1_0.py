import numpy as np
from scipy import integrate

def solve_radiation_ratio():
    """
    Calculates the maximum achievable ratio of bidirectional conical power
    to line intensity for the given rotating charged sphere problem.
    """

    # The intensity function I(theta) is proportional to f(u) where u = cos^2(theta).
    # I(theta) ~ (gamma * sin^2(2*theta) + 1 + cos^2(theta))
    # with gamma = 4 (ratio of E2 to M1 strength for rotation about z-axis).
    # I(theta) ~ 4*(2*sin(theta)*cos(theta))^2 + 1 + cos^2(theta)
    #        = 16*sin^2(theta)*cos^2(theta) + 1 + cos^2(theta)
    # Let u = cos(theta), sin^2(theta) = 1 - u^2.
    # f(u) = 16 * (1 - u**2) * u**2 + 1 + u**2 = -16*u**4 + 17*u**2 + 1.
    # The actual variable for theta is the angle from the z-axis.
    def intensity_function_vs_cos_theta(u):
        # u is cos(theta)
        return -16 * u**4 + 17 * u**2 + 1

    # The line intensity should be the minimum possible intensity.
    # Let's find the minimum of f(u) for u in [-1, 1]. Let x=u^2, x in [0,1].
    # g(x) = -16*x^2 + 17*x + 1.
    # g'(x) = -32x + 17 = 0 -> x = 17/32 (maximum).
    # The minimum is at the boundaries of x in [0,1].
    # At x=0 (theta=pi/2, equator), g(0)=1.
    # At x=1 (theta=0 or pi, poles), g(1)=-16+17+1=2.
    # So, the minimum intensity corresponds to g(0)=1.
    # We choose the line of sight to be in the equatorial plane (theta = pi/2).
    # In our normalized units, I_line = 1.0.
    I_line = 1.0

    # The bidirectional conical power is integrated over two cones with pi/4 half-angle.
    # We align the cones with the z-axis for simplicity and high emission.
    # The integral is over 0 <= theta <= pi/4 and 3pi/4 <= theta <= pi.
    # P_cone = 2 * integral_{2pi} [ integral_{0 to pi/4} I(theta) * sin(theta) d(theta) ] d(phi)
    # Let's define the integrand for the theta integral. u=cos(theta), du=-sin(theta)d(theta)
    integrand = lambda u: intensity_function_vs_cos_theta(u)

    # Integration bounds for the upper cone: theta from 0 to pi/4 -> u from 1 to 1/sqrt(2)
    cos_pi_4 = np.cos(np.pi / 4)
    # We integrate from cos(pi/4) to cos(0)
    integral_upper_cone, _ = integrate.quad(integrand, cos_pi_4, 1)

    # By symmetry, the integral for the lower cone is the same.
    # So, the total integral over theta is 2 * integral_upper_cone
    # P_cone = (Integral over phi=2*pi) * (2 * integral_upper_cone)
    # Don't forget the factor of 2 for the bidirectional cones.
    P_cone = 2 * np.pi * (2 * integral_upper_cone)

    # The final ratio is P_cone / I_line.
    ratio = P_cone / I_line
    
    # We can also calculate this analytically to be sure.
    # integral(-16u^4+17u^2+1)du = -16/5 u^5 + 17/3 u^3 + u
    # Evaluating from 1/sqrt(2) to 1 gives:
    val_at_1 = -16/5 + 17/3 + 1 # = 52/15
    val_at_sqrt2 = (-16/(5*4*np.sqrt(2)) + 17/(3*2*np.sqrt(2)) + 1/np.sqrt(2))
    # = (1/sqrt(2)) * (-4/5 + 17/6 + 1) = (1/sqrt(2)) * ((-24+85+30)/30) = 91/(30*sqrt(2))
    analytical_integral_val = val_at_1 - val_at_sqrt2 # ~ 1.3218
    analytical_p_cone = 4 * np.pi * analytical_integral_val
    
    print(f"The integral from cos(pi/4) to 1 is: {analytical_integral_val}")
    print(f"The cone power (relative units) is 4 * pi * {analytical_integral_val} = {analytical_p_cone}")
    print(f"The minimum line intensity (relative units) is: {I_line}")
    print(f"The maximum achievable ratio is P_cone / I_line = {analytical_p_cone} / {I_line}")
    print("Final equation: P_cone / I_line = (4 * pi * (-16/5 * (1)^5 + 17/3 * (1)^3 + 1 - (-16/5 * (1/sqrt(2))^5 + 17/3 * (1/sqrt(2))^3 + 1/sqrt(2)))) / 1")
    # For a clean final answer
    print("\nFinal calculation:")
    # Using more explicit numbers from the formula
    P_cone_num = 4 * np.pi * (52/15 - 91/(30*np.sqrt(2)))
    print(f"Ratio = 4 * pi * (52/15 - 91/(30 * sqrt(2))) = {P_cone_num}")


solve_radiation_ratio()