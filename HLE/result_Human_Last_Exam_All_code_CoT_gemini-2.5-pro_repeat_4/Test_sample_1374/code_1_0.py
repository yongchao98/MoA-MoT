import math

def calculate_max_distance():
    """
    Calculates the furthest distance from point A for a playdough shape
    optimized to produce the strongest gravitational field.
    """
    # The total volume of the playdough in cubic meters.
    volume = 1.0

    # The formula for the maximum distance 'r_max' is derived from optimizing
    # the gravitational field and is given by:
    # r_max = (15 * V / (4 * pi))^(1/3)

    # Breaking down the calculation for clarity as requested.
    # Numerator of the base of the exponent
    numerator = 15 * volume
    # Denominator of the base of the exponent
    denominator = 4 * math.pi
    # The base of the exponent
    base = numerator / denominator
    # The exponent
    exponent = 1.0 / 3.0

    # Calculate the final result
    r_max = base ** exponent

    # Output the equation with the numbers plugged in
    print("The equation for the furthest distance is derived as r_max = ( (15 * V) / (4 * pi) )^(1/3)")
    print(f"Plugging in the numbers: r_max = ( ({numerator}) / ({denominator:.4f}) )^({exponent:.3f})")

    # Output the final answer
    print(f"\nThe furthest point on the surface of the playdough is {r_max:.6f} meters from point A.")

calculate_max_distance()