import math

def calculate_slide_angle(L, l, mu):
    """
    Calculates the angle theta at which a tilting rod begins to slide.

    Args:
        L (float): The total length of the rod.
        l (float): The distance from the pivot to the center of mass.
        mu (float): The coefficient of static friction.

    Returns:
        float: The angle theta in degrees.
    """
    # Check for physical constraints
    if l >= L / 2:
        print("Error: l must be less than L/2 for the setup to be physical.")
        return None
    if L <= 0 or l < 0 or mu < 0:
        print("Error: Lengths and friction coefficient must be positive.")
        return None

    # The derived relationship is tan(theta) = mu * (L^2 + 12*l^2) / (L^2 + 36*l^2)
    
    # Calculate terms
    L_sq = L**2
    l_sq = l**2
    
    numerator = L_sq + 12 * l_sq
    denominator = L_sq + 36 * l_sq
    
    # Calculate tan(theta)
    tan_theta = mu * (numerator / denominator)
    
    # Calculate theta in radians and then convert to degrees
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)
    
    # Print the equation with the numbers plugged in
    print("The equation for the slide angle is: tan(theta) = mu * (L^2 + 12*l^2) / (L^2 + 36*l^2)")
    print("\nPlugging in the values L={}, l={}, mu={}:".format(L, l, mu))
    print("tan(theta) = {} * (({}^2 + 12*{}^2) / ({}^2 + 36*{}^2))".format(mu, L, l, L, l))
    print("tan(theta) = {} * (({} + {}) / ({} + {}))".format(mu, round(L_sq, 4), round(12 * l_sq, 4), round(L_sq, 4), round(36 * l_sq, 4)))
    print("tan(theta) = {} * ({} / {})".format(mu, round(numerator, 4), round(denominator, 4)))
    print("tan(theta) = {}".format(round(tan_theta, 4)))
    print("\nThis gives the final angle:")
    print("theta = arctan({}) = {} degrees".format(round(tan_theta, 4), round(theta_deg, 2)))
    
    return theta_deg

# --- Example Usage ---
# You can change these values to see how the angle changes.
rod_length_L = 1.0  # meters
cm_distance_l = 0.1 # meters (distance from pivot to CM)
friction_coeff_mu = 0.5 # dimensionless

# Execute the calculation
calculate_slide_angle(rod_length_L, cm_distance_l, friction_coeff_mu)
