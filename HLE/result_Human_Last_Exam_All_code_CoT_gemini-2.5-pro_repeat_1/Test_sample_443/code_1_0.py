import numpy as np

def solve():
    """
    This function demonstrates the lower bound for k by constructing a
    polynomial for which the number of required balls grows as D^2.
    """
    # The degree D of the polynomial P.
    D = 10

    # The angle constraint is > 1/10 radians.
    angle_rad = 0.1
    c = np.sin(angle_rad)

    # We constructed a polynomial P whose zero set Z(P,T) has a
    # vertical extent of 2*L, where L is proportional to D^2.
    # The constant of proportionality is derived from the angle constraint.
    # We choose a valid L = factor * D^2.
    # The upper limit for L is (sqrt(1-c^2)/(4*c)) * D^2. We take half of that.
    proportionality_constant = (np.sqrt(1 - c**2) / (4 * c)) * 0.5
    
    L = proportionality_constant * D**2

    # The vertical extent of the constructed surface inside the cylinder is 2*L.
    vertical_extent = 2 * L

    # The number of unit balls (diameter 2) needed to cover this vertical extent
    # is at least the extent divided by the ball diameter.
    num_balls = vertical_extent / 2.0

    print("Demonstration for the lower bound k >= 2:")
    print(f"For a polynomial of degree D = {D}, we construct an example surface.")
    print(f"The number of unit balls needed to cover this surface is at least {num_balls:.2f}.")
    print("\nThis number is proportional to D^k.")
    print("From our construction, the equation for the number of balls is:")
    
    k = 2
    # The prompt asks to output each number in the final equation.
    # Equation: Num_Balls = proportionality_constant * D ^ k
    print(f"Num_Balls = {proportionality_constant:.4f} * {D}^{k}")

solve()
