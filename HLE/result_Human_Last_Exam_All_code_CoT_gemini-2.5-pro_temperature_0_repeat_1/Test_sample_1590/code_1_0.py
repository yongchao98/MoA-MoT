import math

def calculate_sliding_angle(L, l, mu):
    """
    Calculates the angle theta at which the rod begins to slide.

    Args:
        L (float): The total length of the rod.
        l (float): The distance from the pivot (table edge) to the rod's center of mass.
        mu (float): The coefficient of static friction between the rod and the table edge.

    Returns:
        float: The angle theta in degrees.
    """
    # The expression for tan(theta) is derived from the condition f = mu * N,
    # including the effects of rotational motion.
    # tan(theta) = (mu*L^2/12 + l^2*(1+mu)) / (L^2/12 + l^2*(1+2*mu))

    # Calculate the terms in the expression
    L_sq_div_12 = L**2 / 12
    l_sq = l**2

    numerator = mu * L_sq_div_12 + l_sq * (1 + mu)
    denominator = L_sq_div_12 + l_sq * (1 + 2 * mu)

    # Avoid division by zero if the denominator is zero
    if denominator == 0:
        return float('inf') # Or handle as an error

    tan_theta = numerator / denominator
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)
    
    return theta_deg, numerator, denominator, tan_theta

# --- Example Parameters ---
# A rod of length L=1.0m lies on a table.
# The length hanging off is l + L/2. Let's assume l = 0.1m.
# So, 0.1 + 1.0/2 = 0.6m hangs off, and 0.4m is on the table.
# The pivot is the edge, and the center of mass is at 0.5m from the end,
# so the distance 'l' from the pivot to the CM is 0.5m - 0.4m = 0.1m.
L_rod = 1.0  # meters
l_cm_dist = 0.1 # meters
mu_friction = 0.5 # coefficient of static friction

# --- Calculation and Output ---
angle_deg, num, den, tan_val = calculate_sliding_angle(L_rod, l_cm_dist, mu_friction)

print("The general expression for the angle theta is:")
print("theta = arctan( (mu*L^2/12 + l^2*(1+mu)) / (L^2/12 + l^2*(1+2*mu)) )")
print("\n--- Calculation for Specific Values ---")
print(f"Rod Length (L): {L_rod} m")
print(f"CM distance from pivot (l): {l_cm_dist} m")
print(f"Coefficient of Friction (mu): {mu_friction}")

print("\nPlugging these numbers into the equation for tan(theta):")
# Using f-strings to display the equation with the numbers
print(f"tan(theta) = ({mu_friction} * {L_rod**2:.2f}/12 + {l_cm_dist**2:.2f} * (1 + {mu_friction})) / ({L_rod**2:.2f}/12 + {l_cm_dist**2:.2f} * (1 + 2 * {mu_friction}))")
print(f"tan(theta) = ({num:.4f}) / ({den:.4f})")
print(f"tan(theta) = {tan_val:.4f}")

print("\nSolving for the angle:")
print(f"theta = arctan({tan_val:.4f})")
print(f"theta = {angle_deg:.2f} degrees")
