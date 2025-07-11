import numpy as np

def solve_retarded_time(v, c, L, d, is_front):
    """
    Solves the quadratic equation for the retarded time t_e.
    (c^2 - v^2) * t^2 - (2*v*x0) * t - (x0^2 + d^2) = 0
    t_obs is assumed to be 0. We need the negative root t_e < 0.
    """
    x0 = L / 2 if is_front else -L / 2
    
    a = c**2 - v**2
    b = -2 * v * x0
    const_c = -(x0**2 + d**2)
    
    # Quadratic formula: t = [-b +/- sqrt(b^2 - 4ac)] / 2a
    discriminant = b**2 - 4 * a * const_c
    if discriminant < 0:
        return None # No real solution
    
    # We need the negative time solution since t_e < t_obs=0
    t_e = (-b - np.sqrt(discriminant)) / (2 * a)
    return t_e

def calculate_shift():
    """
    Calculates and prints the shift in the apparent center of gravity.
    """
    # Parameters
    c = 299792458.0  # Speed of light in m/s
    v = 0.5 * c      # Velocity of mass 2 (50% of c)
    d = 1.0e9        # Distance of closest approach in meters
    L = 1.0e8        # Length of mass 2 in meters

    # --- Step 1: Calculate retarded times and positions ---
    # For the front of the object
    t_ef = solve_retarded_time(v, c, L, d, is_front=True)
    pos_f = np.array([v * t_ef + L / 2, d])

    # For the back of the object
    t_eb = solve_retarded_time(v, c, L, d, is_front=False)
    pos_b = np.array([v * t_eb - L / 2, d])
    
    # --- Step 2: Calculate retarded distances ---
    R_f = np.linalg.norm(pos_f)
    R_b = np.linalg.norm(pos_b)

    # --- Step 3: Apply Assumption C and calculate forces ---
    # Assumption C: Field strength (magnitude) is proportional to 1/T_app = c/R
    # We can use a proportionality constant of 1 for simplicity.
    mag_f = 1.0 / R_f
    mag_b = 1.0 / R_b

    # Force vector = magnitude * unit_direction_vector
    # Unit direction from retarded position P to observer at origin is -P/R
    force_f = mag_f * (-pos_f / R_f)
    force_b = mag_b * (-pos_b / R_b)
    
    total_force = force_f + force_b

    # --- Step 4: Determine apparent positions and compare ---
    # The "default" apparent position would be the geometric center of the retarded positions
    retarded_center_pos = (pos_f + pos_b) / 2.0
    
    # The new apparent position is in the direction opposite to the total force vector.
    # We only need the direction vector for comparison.
    apparent_source_dir = -total_force / np.linalg.norm(total_force)
    
    # For a meaningful comparison, we project this direction onto the y=d plane
    # The apparent x position is where this line of force crosses the y=d line
    apparent_x_pos = d * (apparent_source_dir[0] / apparent_source_dir[1])
    
    print("Demonstrating the effect of Assumption C:")
    print(f"Object m2 is moving in the +x direction with velocity v = {v/c:.2f}c")
    print("-" * 50)
    print("Retarded position of the front: ({:.3e}, {:.3e}) m".format(pos_f[0], pos_f[1]))
    print("Retarded position of the back:  ({:.3e}, {:.3e}) m".format(pos_b[0], pos_b[1]))
    print("-" * 50)
    
    # Calculate x-position of the retarded geometric center
    retarded_center_x = retarded_center_pos[0]
    print(f"X-position of the RETARDED GEOMETRIC CENTER: {retarded_center_x:,.2f} m")
    
    # Calculate the apparent x-position based on the force calculation
    print(f"X-position of the APPARENT CENTER OF GRAVITY: {apparent_x_pos:,.2f} m")
    print("-" * 50)

    shift = apparent_x_pos - retarded_center_x
    if shift > 0:
        print(f"The apparent center of gravity is shifted by {shift:,.2f} m.")
        print("This shift is in the direction of motion (+x), as predicted.")
    else:
        print("The apparent center of gravity is shifted opposite to the direction of motion.")

calculate_shift()
<<<C>>>