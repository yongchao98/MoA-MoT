import math

def solve_circle_tangency():
    """
    This function solves for r^2 based on the geometric constraints.
    The problem simplifies to solving the quadratic equation r^2 - 4r + 4 = 0.
    """
    
    # Coefficients of the quadratic equation ar^2 + br + c = 0
    a = 1
    b = -4
    c = 4

    # The equation derived from the geometric setup is (r+2)^2 = 2*r^2 + 8,
    # which simplifies to r^2 - 4r + 4 = 0.
    
    # We solve this quadratic equation for r.
    # The discriminant is d = b**2 - 4*a*c
    discriminant = b**2 - 4*a*c

    # Since the discriminant is 0, there is exactly one real solution for r.
    r = -b / (2 * a)

    # Calculate the value of r^2
    r_squared = r**2

    print("The geometric analysis leads to the following quadratic equation for the radius r:")
    print(f"({a})r^2 + ({b})r + ({c}) = 0")
    print(f"Solving this equation gives a unique positive radius: r = {int(r)}")
    print(f"The value of r^2 is {int(r)}^2, which is:")
    print(int(r_squared))

solve_circle_tangency()