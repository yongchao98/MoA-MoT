import numpy as np

def calculate_angular_distance():
    """
    Calculates the angular distance between two stars based on precession data.

    The problem provides the following information:
    - Precession period T = 26000 years.
    - Earth's axial tilt epsilon = 23.5 degrees.
    - Star A was last on the celestial equator 3000 years ago.
    - Star B will first be on the celestial equator in 10000 years.
    - After a certain time, the stars will swap equatorial coordinates.
    - Both stars are currently on the same side of the celestial equator.

    The solution proceeds as follows:
    1. The coordinate swap implies the stars' ecliptic coordinates are related. Analysis shows
       two possibilities:
       a) beta_A = beta_B, lambda_A = lambda_B + 180 deg
       b) beta_A = -beta_B, lambda_A = lambda_B
       Further analysis of the equator crossing directions shows that only model (b) is possible.

    2. For model (b), the angular distance `theta` is given by cos(theta) = cos(2*beta).
       So, theta = 2 * |beta|.

    3. The declination of the stars can be modeled as:
       sin(delta_A) = cos(epsilon)*sin(beta) + sin(epsilon)*cos(beta)*cos(omega*t + phi)
       sin(delta_B) = -cos(epsilon)*sin(beta) + sin(epsilon)*cos(beta)*cos(omega*t + phi)
       where beta = beta_A.

    4. The equator crossing times are t_A = -3000 and t_B = 10000.
       The time difference is 13000 years, which is exactly T/2. This provides a consistency
       check. The specific values suggest that t=0 is a time of extremum for the cosine term,
       meaning the phase phi is 0 or pi. We assume phi=0.

    5. With phi=0, we can solve for beta:
       tan(beta) / tan(epsilon) = cos(omega * t_B)
       tan(beta) / tan(epsilon) = -cos(omega * t_A)
       The values are consistent: cos(omega*10000) = -cos(omega*3000).

    6. We calculate beta and then the angular distance theta.
    """
    # Given parameters
    T = 26000.0  # years
    epsilon_deg = 23.5
    epsilon_rad = np.deg2rad(epsilon_deg)
    
    t_A = -3000.0  # years
    t_B = 10000.0 # years
    
    # Angular frequency of precession
    omega = 2 * np.pi / T  # rad/year
    
    # From the problem's constraints, we deduce that the phase of the precessional
    # cycle is zero at the current time (t=0). This gives a specific value for
    # the ratio tan(beta)/tan(epsilon).
    # tan(beta) / tan(epsilon) = cos(omega * t_B)
    cos_val = np.cos(omega * t_B)
    
    # Calculate tan(beta)
    tan_beta = np.tan(epsilon_rad) * cos_val
    
    # Calculate beta in radians
    beta_rad = np.arctan(tan_beta)
    
    # The angular distance theta for our derived model (beta_A = -beta_B, lambda_A = lambda_B)
    # is 2 * |beta|.
    theta_rad = 2 * np.abs(beta_rad)
    
    # Convert final answer to degrees
    theta_deg = np.rad2deg(theta_rad)
    
    # Output the logic and the final result
    print("Step 1: Determine the relationship between the stars' ecliptic coordinates.")
    print("The coordinate swap condition and equator crossing data imply that the stars have the same ecliptic longitude and opposite ecliptic latitudes (beta_B = -beta_A).")
    print("\nStep 2: Relate the angular distance to the ecliptic latitude.")
    print("The angular distance theta between them is given by theta = 2 * |beta_A|.")
    print("\nStep 3: Use the equator crossing times to find the value of beta_A.")
    print("The declination model is sin(delta) = C + A*cos(omega*t).")
    print(f"The given times t_A = {int(t_A)} years and t_B = {int(t_B)} years are used.")
    print("This leads to the equation: tan(beta_A) = tan(epsilon) * cos(omega * t_B).")
    print(f"tan(beta_A) = tan({epsilon_deg}) * cos(2*pi/{int(T)} * {int(t_B)})")
    
    # Calculation part for printing
    tan_epsilon_val = np.tan(epsilon_rad)
    beta_deg = np.rad2deg(beta_rad)
    
    print(f"tan(beta_A) = {tan_epsilon_val:.4f} * {cos_val:.4f} = {tan_beta:.4f}")
    print(f"beta_A = arctan({tan_beta:.4f}) = {beta_deg:.2f} degrees.")
    
    print("\nStep 4: Calculate the final angular distance.")
    print(f"theta = 2 * |{beta_deg:.2f}| = {theta_deg:.2f} degrees.")
    
    # Round to one decimal place for the final answer
    final_answer = round(theta_deg, 1)
    print(f"\nThe angular distance between the two stars is {final_answer} degrees.")

calculate_angular_distance()