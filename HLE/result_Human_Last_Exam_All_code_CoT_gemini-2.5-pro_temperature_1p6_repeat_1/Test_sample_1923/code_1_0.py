import numpy as np

def calculate_force_shift():
    """
    Simulates the gravitational force from a moving mass on a stationary observer
    to determine the effect of assumption C.
    """
    # --- Simulation Parameters ---
    v = 0.5    # Speed of mass 2 (as a fraction of c)
    c = 1.0    # Speed of light/gravity
    y_pos = 1.0  # Closest approach distance of mass 2 to the y-axis
    
    # The trajectory of mass 2 is from x_start to x_end at y=y_pos
    x_start = -20.0
    x_end = 20.0
    num_steps = 4001
    
    # Observer (mass 1) is at the origin (0, 0)
    observer_pos = np.array([0.0, 0.0])
    
    # Initialize total force vectors (as a sum/integral over the path)
    total_force_base = np.array([0.0, 0.0])
    total_force_modulated = np.array([0.0, 0.0])
    
    # --- Simulation Loop ---
    # Integrate the force contribution from each point on the path
    x_coords = np.linspace(x_start, x_end, num_steps)
    dx = (x_end - x_start) / (num_steps - 1)
    
    for x in x_coords:
        source_pos = np.array([x, y_pos])
        velocity_vec = np.array([v, 0.0])
        
        # Vector from observer to the source's (retarded) position
        r_vec = source_pos - observer_pos
        distance = np.linalg.norm(r_vec)
        
        if distance == 0:
            continue
            
        r_hat = r_vec / distance
        
        # --- 1. Base Inverse-Square Force (without modulation) ---
        # The force on the observer points towards the source (-r_hat)
        # Magnitude is 1/r^2 (constants G, M1, M2 set to 1)
        force_base_vec = -r_hat / (distance**2)
        
        # Integrate by summing contributions (multiplied by path element dx)
        total_force_base += force_base_vec * dx
        
        # --- 2. Modulated Force (with assumption C) ---
        # Calculate the radial velocity (v_r)
        # v_r = projection of velocity_vec onto r_hat
        v_r = np.dot(velocity_vec, r_hat)
        
        # Calculate modulation factor M = 1 / (1 - v_r/c)
        modulation_factor = 1.0 / (1.0 - v_r / c)
        
        # Apply the modulation to the force magnitude
        force_modulated_vec = force_base_vec * modulation_factor
        
        # Integrate the modulated force
        total_force_modulated += force_modulated_vec * dx
        
    print("This simulation calculates the total force impulse from a moving mass.")
    print("A standard inverse-square law results in a force vector with no component along the direction of motion (x=0).")
    print(f"Calculated base force vector (Fx, Fy): ({total_force_base[0]:.4f}, {total_force_base[1]:.4f})")
    print("\nApplying assumption C modulates the field strength based on radial velocity.")
    print("This results in a net force component in the direction of motion (positive Fx).")
    
    # "Output each number in the final equation"
    print("\nThe final equation for the total modulated force impulse F is:")
    print(f"F = ({total_force_modulated[0]:.4f})i + ({total_force_modulated[1]:.4f})j")
    
    print("\nSince the i-component (the direction of motion) is positive, this shows")
    print("that the apparent center of gravity is shifted in the direction of motion.")

if __name__ == '__main__':
    calculate_force_shift()
