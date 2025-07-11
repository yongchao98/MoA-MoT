import math

def find_sliding_angle(mu):
    """
    Calculates the angle at which the rod begins to slide based on the 
    coefficient of static friction (mu).

    Args:
        mu (float): The coefficient of static friction.

    Returns:
        float: The sliding angle in degrees, or None if input is invalid.
    """
    if mu < 0:
        print("Error: The coefficient of friction cannot be negative.")
        return None
    
    # The derivation shows that the condition for sliding is tan(theta) = mu.
    # From this, we can find theta.
    # The final equation is theta = arctan(mu).
    
    # Calculate theta in radians using math.atan()
    theta_radians = math.atan(mu)
    
    # Convert theta to degrees
    theta_degrees = math.degrees(theta_radians)
    
    print("--- Derivation Summary ---")
    print("1. The normal force on the rod is N = M*g*cos(theta).")
    print("2. The force pulling the rod down its slope is F_parallel = M*g*sin(theta).")
    print("3. The maximum static friction force is f_max = mu * N = mu*M*g*cos(theta).")
    print("4. Sliding begins when F_parallel = f_max.")
    print("   M*g*sin(theta) = mu*M*g*cos(theta)")
    print("5. Simplifying gives the final equation: tan(theta) = mu")
    print("--------------------------\n")
    
    print(f"For a coefficient of friction mu = {mu}:")
    print(f"The equation is tan(theta) = {mu}")
    print(f"The angle theta at which sliding begins is arctan({mu}).")
    print(f"Result: theta ≈ {theta_degrees:.2f} degrees ({theta_radians:.4f} radians).")
    
    return theta_degrees

if __name__ == '__main__':
    # You can change this value to any non-negative number to see the result.
    coefficient_of_friction = 0.75
    find_sliding_angle(coefficient_of_friction)
