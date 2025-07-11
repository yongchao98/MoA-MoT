import math

def simulate_gravity_shift():
    """
    Simulates the gravitational force from a moving rod on a stationary observer
    to demonstrate the principle described in the problem.
    """
    # --- Parameters ---
    G = 6.674e-11  # Gravitational constant (not needed for direction)
    M2 = 1.0         # Total mass of the moving rod (mass 2)
    L = 2.0          # Length of the rod
    Y = 5.0          # Perpendicular distance to the rod's path
    v = 0.5          # Velocity of the rod as a fraction of c
    c = 1.0          # Speed of light/gravity
    N = 200          # Number of segments to model the rod

    # Mass of each segment
    dm = M2 / N
    # Position of observer (mass 1)
    m1_pos = (0, 0)

    # Initialize total force vector
    total_force_x = 0.0
    total_force_y = 0.0

    # Integrate force over the length of the rod
    for i in range(N):
        # Position of the mass segment 'dm' on the rod
        # Rod extends from -L/2 to L/2
        x = -L/2 + (i + 0.5) * (L / N)
        dm_pos = (x, Y)

        # --- Calculate force from this segment ---
        
        # Vector from observer to source segment
        r_vec_x = dm_pos[0] - m1_pos[0]
        r_vec_y = dm_pos[1] - m1_pos[1]
        
        # Distance to the segment
        R = math.sqrt(r_vec_x**2 + r_vec_y**2)

        # Unit vector for force direction (attractive: from observer to source)
        u_vec_x = r_vec_x / R
        u_vec_y = r_vec_y / R

        # --- Apply Assumption C ---
        # "Field strength varies inversely with apparent propagation time"
        # We model this as a relativistic beaming effect that strengthens the field
        # from the leading part of the object.
        # This requires a hypothetical force law S ~ 1/(1 - v*cos(theta)/c)
        # where theta is the angle between v and the vector to the source.
        # cos(theta) = x / R
        beta = v / c
        # The modulation factor for field strength 'S'
        # This specific form makes the leading edge (x>0) stronger
        strength_factor = 1.0 / (1.0 - beta * (x / R))
        
        # Magnitude of the force from the segment, modulated by the strength factor
        # dF = S * (G * dm * m1 / R^2), m1=1
        dF_mag = strength_factor * (G * dm / R**2)

        # Force vector from the segment
        dF_x = dF_mag * u_vec_x
        dF_y = dF_mag * u_vec_y

        # Add to total force
        total_force_x += dF_x
        total_force_y += dF_y

    # --- Output the results ---
    print("Demonstration of principle from Option C:")
    print("-" * 40)
    print("A simulation was run with the following parameters:")
    print(f"  Velocity (v/c): {v}")
    print(f"  Rod Length (L): {L} m")
    print(f"  Observer Distance (Y): {Y} m")
    print("\nBased on the assumption that field strength is modulated by motion,")
    print("we calculate the total gravitational force vector on the observer.")
    print("The final equation for the total force F = (Fx, Fy) is a sum over rod segments:")
    print("F = Σ [ S(x) * (G*dm/R²) * (x/R, Y/R) ]")
    print(f"where the strength factor S(x) = 1 / (1 - ({v:.2f}) * (x/R))")
    
    print("\n--- Simulation Results ---")
    print(f"Calculated Total Force (Fx): {total_force_x / G:.4f} * G")
    print(f"Calculated Total Force (Fy): {total_force_y / G:.4f} * G")
    print("-" * 40)
    print("\nConclusion:")
    print("The x-component of the force (Fx) is positive.")
    print("Since the rod is moving in the positive x-direction, this force component")
    print("pulls the observer towards the leading side of the rod.")
    print("This means the center of gravity of the moving rod appears shifted")
    print("in the direction of its motion.")

if __name__ == '__main__':
    simulate_gravity_shift()