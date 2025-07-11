import math

def calculate_max_distance():
    """
    Calculates the maximum distance from point A to the surface of a 1 cubic meter
    object shaped to create the strongest possible gravitational field at A.
    """
    
    # The volume of the playdough in cubic meters.
    volume = 1.0
    
    # The formula for the furthest point (r_max) is derived from optimizing the 
    # gravitational field integral.
    # r_max = (15 * V / (4 * pi))^(1/3)
    
    # Define the numbers in the final equation
    numerator_const = 15.0
    denominator_const = 4.0
    pi_val = math.pi
    exponent = 1.0 / 3.0
    
    # Calculate the result
    numerator = numerator_const * volume
    denominator = denominator_const * pi_val
    result = (numerator / denominator) ** exponent

    print("To find the furthest point on the surface, we use the equation:")
    print(f"r_max = ({int(numerator_const)} * V / ({int(denominator_const)} * pi)) ^ (1/3)")
    print("\nPlugging in the values:")
    print(f"V = {volume} m^3")
    print(f"pi = {pi_val}")
    print("\nFinal Calculation:")
    print(f"r_max = ({numerator} / ({denominator_const} * {pi_val})) ^ ({exponent})")
    print(f"r_max = {result} meters")

calculate_max_distance()