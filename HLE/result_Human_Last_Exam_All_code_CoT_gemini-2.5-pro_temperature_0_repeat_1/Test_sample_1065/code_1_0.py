import math

def calculate_string_center_of_mass():
    """
    Calculates the coordinates of the center of mass of the candy string.

    The problem is solved by first finding the equilibrium position of the string on
    the one-fourth spherical pumpkin. The string consists of two parts:
    1. An arc on the sphere's surface.
    2. A vertical segment hanging freely.

    The calculation proceeds as follows:
    - The lift-off angle (theta_0) where the string leaves the sphere is determined
      by force balance, resulting in cos(theta_0) = 1/2, so theta_0 = pi/3.
    - The length of the hanging part is determined to be R/2.
    - We set the radius R=1 for a numerical solution.
    - The center of mass of the arc and the hanging line are calculated separately.
    - The total center of mass is found using the composite body formula.
    """
    # Set radius R to 1 for numerical calculation.
    R = 1.0
    
    # Lift-off angle in radians
    theta_0 = math.pi / 3.0
    
    # Part 1: The arc on the sphere
    # Length of the arc (proportional to its mass)
    L1 = R * theta_0
    # Center of mass of the arc
    # x_cm1 = R * (1 - cos(theta_0)) / theta_0
    x1 = R * (1 - math.cos(theta_0)) / theta_0
    # z_cm1 = R * sin(theta_0) / theta_0
    z1 = R * math.sin(theta_0) / theta_0
    
    # Part 2: The hanging vertical line
    # Length of the hanging line (proportional to its mass)
    L2 = R / 2.0
    # Center of mass of the hanging line
    # The line hangs from x_p = R * sin(theta_0)
    x2 = R * math.sin(theta_0)
    # The line hangs from z_p = R * cos(theta_0) = R/2 down to z=0.
    # The center of mass is at the midpoint, z = (R/2)/2 = R/4.
    z2 = R / 4.0
    
    # Total length (proportional to total mass)
    L_total = L1 + L2
    
    # Calculate the combined center of mass
    # X_cm = (L1*x1 + L2*x2) / (L1 + L2)
    X_cm = (L1 * x1 + L2 * x2) / L_total
    # Z_cm = (L1*z1 + L2*z2) / (L1 + L2)
    Z_cm = (L1 * z1 + L2 * z2) / L_total
    
    # Print the components of the calculation as requested
    print("--- Calculation Details (assuming R=1) ---")
    print(f"Part 1 (Arc):")
    print(f"  Length (Mass) M1 = {L1}")
    print(f"  Center of Mass (x1, z1) = ({x1}, {z1})")
    print(f"Part 2 (Hanging Line):")
    print(f"  Length (Mass) M2 = {L2}")
    print(f"  Center of Mass (x2, z2) = ({x2}, {z2})")
    print(f"Total Mass M_total = {L_total}")
    print("\n--- Final Coordinates ---")
    print("The horizontal and vertical coordinates of the center of mass are:")
    # The final output is printed as a comma-separated string
    print(f"{X_cm},{Z_cm}")

calculate_string_center_of_mass()