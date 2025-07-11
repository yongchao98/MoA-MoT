import math

def calculate_max_distance():
    """
    Calculates the furthest distance from point A to the surface of a playdough shape
    optimized to produce the strongest gravitational field at A.

    The formula for this distance is r_max = (15 / (4 * pi))^(1/3).
    """

    # Numbers in the final equation
    numerator = 15
    denominator_const = 4
    power_num = 1
    power_den = 3

    print("The final equation for the maximum distance (r_max) is:")
    print(f"r_max = ({numerator} / ({denominator_const} * pi)) ^ ({power_num}/{power_den})")
    print("\nHere are the values used in the calculation:")
    print(f"Numerator: {numerator}")
    print(f"Constant in denominator: {denominator_const}")
    print(f"Value of pi: {math.pi}")
    print(f"Exponent: {power_num}/{power_den}")

    # Perform the calculation
    base = numerator / (denominator_const * math.pi)
    exponent = power_num / power_den
    result = base ** exponent

    print(f"\nThe calculated furthest distance is:")
    print(result)

if __name__ == "__main__":
    calculate_max_distance()
