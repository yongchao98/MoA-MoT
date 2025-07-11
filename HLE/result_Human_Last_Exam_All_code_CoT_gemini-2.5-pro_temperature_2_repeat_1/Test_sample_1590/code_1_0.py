import math

def calculate_sliding_angle(L, l, mu):
    """
    Calculates and displays the angle theta at which a rod begins to slide.

    The function first shows the derived symbolic formula, then substitutes
    the provided numerical values to calculate the angle, showing the intermediate steps.

    Args:
        L (float): The total length of the rod.
        l (float): The distance from the rod's center of mass to the pivot point (table edge).
        mu (float): The coefficient of static friction.
    """
    # Verify inputs are physically reasonable
    if not (L > 0 and l >= 0 and mu >= 0):
        print("Error: Physical parameters (L, l, mu) must be non-negative, and L must be positive.")
        return
    if (l > L / 2):
        print(f"Error: The distance l ({l}) cannot be greater than half the rod length ({L/2}).")
        return

    print("The derived expression for the angle theta where the rod begins to slide is:")
    print("tan(theta) = mu * L^2 / (L^2 + 36 * l^2)\n")

    print("--- Calculation with provided values ---")
    print(f"Length of the rod (L) = {L}")
    print(f"Distance of CM from pivot (l) = {l}")
    print(f"Coefficient of friction (mu) = {mu}\n")
    
    # Calculate terms for the equation
    L_squared = L**2
    l_squared = l**2
    
    numerator = mu * L_squared
    denominator = L_squared + 36 * l_squared
    
    # Display the final equation with numerical values, showing each component
    print("Substituting these values into the equation:")
    # Printing with f-string formatting to show the equation construction
    print(f"tan(theta) = {mu} * ({L})^2 / (({L})^2 + 36 * ({l})^2)")
    print(f"tan(theta) = {mu} * {L_squared} / ({L_squared} + 36 * {l_squared})")
    print(f"tan(theta) = {numerator} / ({L_squared} + {36 * l_squared})")
    print(f"tan(theta) = {numerator} / {denominator}\n")
    
    # Calculate tan(theta)
    if denominator == 0:
        # This case is physically unlikely unless L=l=0
        print("Error: Denominator is zero, unable to calculate angle.")
        return
        
    tan_theta = numerator / denominator
    print(f"Resulting tan(theta) = {tan_theta:.4f}\n")
    
    # Calculate the angle in degrees
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)
    
    print("--- Final Answer ---")
    print(f"The angle at which the rod begins to slide is theta = atan({tan_theta:.4f})")
    print(f"theta ≈ {theta_deg:.2f} degrees")


if __name__ == '__main__':
    # Example Usage:
    # Define some plausible physical parameters for the rod problem.
    # L: total length of the rod in meters.
    # l: distance of the center of mass from the pivot in meters.
    #    (This means a length of L/2 + l hangs off the edge).
    # mu: coefficient of static friction (dimensionless).
    
    rod_length = 1.0  # meters
    cm_distance_from_edge = 0.1 # meters
    friction_coefficient = 0.5

    calculate_sliding_angle(rod_length, cm_distance_from_edge, friction_coefficient)