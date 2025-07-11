import math

def solve_gravity_shape_problem():
    """
    Calculates the furthest point on a 1m^3 object shaped to maximize gravity.
    """
    # The problem asks for the furthest point on a shape of volume 1m^3
    # that maximizes the gravitational field at a point A.
    # The optimal shape is a solid of revolution with surface r = R * sqrt(cos(theta)).
    # Its volume V is given by the formula V = (4 * pi * R^3) / 15.

    # We are given V = 1. We solve for R, which is the maximum distance.
    # 1 = (4 * pi * R^3) / 15
    # R^3 = 15 / (4 * pi)
    # R = (15 / (4 * pi))^(1/3)

    # Define the numbers in the equation
    volume = 1.0
    numerator = 15
    denominator_const = 4
    pi_val = math.pi
    
    # Perform the calculation
    r_cubed = numerator / (denominator_const * pi_val)
    r_max = r_cubed**(1/3)

    print("The problem is to find the maximum distance from point A to the surface of a 1 m^3 body shaped to maximize the gravitational field at A.")
    print("The formula for the maximum distance R is derived from the volume V of the optimal shape:")
    print("V = (4 * pi * R^3) / 15")
    print("\nSetting V = 1 and solving for R, we get:")
    print(f"R = ({numerator} / ({denominator_const} * pi))^(1/3)")
    
    # Output the final numerical result
    print(f"\nThe furthest distance from A is approximately {r_max:.10f} meters.")

    # Final answer in the required format
    print(f"<<<{r_max}>>>")

solve_gravity_shape_problem()