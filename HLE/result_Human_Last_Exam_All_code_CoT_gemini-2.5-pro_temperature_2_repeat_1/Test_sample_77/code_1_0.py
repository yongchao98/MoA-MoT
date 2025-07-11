import time

def simulate_fluid_particle(start_x, start_y):
    """
    Simulates a single fluid particle falling in the described scene.
    Prints the outcome based on its starting position.
    """
    # --- Simulation Setup ---
    # These values represent the physical setup you described
    INFLOW_POSITION = {'x': start_x, 'y': start_y}
    OBSTACLE_Y_LEVEL = 2.0
    OBSTACLE_X_MIN = -5.0
    OBSTACLE_X_MAX = 5.0
    DOMAIN_FLOOR_Y = 0.0
    GRAVITY = 0.5  # Represents downward acceleration per time step

    # --- Particle Initialization ---
    particle = INFLOW_POSITION.copy()
    print(f"\n--- Simulating new particle starting at (x={particle['x']}, y={particle['y']}) ---")
    print("Particle emitted from inflow object.")
    
    # --- Simulation Loop ---
    is_falling = True
    while is_falling:
        # Apply gravity
        particle['y'] -= GRAVITY

        # Check for collision with the obstacle
        # Condition: Particle is at or below the obstacle's height AND within its x-bounds
        if particle['y'] <= OBSTACLE_Y_LEVEL and OBSTACLE_X_MIN <= particle['x'] <= OBSTACLE_X_MAX:
            particle['y'] = OBSTACLE_Y_LEVEL
            print(f"Particle hit the obstacle plane at y={OBSTACLE_Y_LEVEL:.1f} and came to rest.")
            is_falling = False # Stop simulation for this particle
            
        # Check for collision with the domain floor
        elif particle['y'] <= DOMAIN_FLOOR_Y:
            particle['y'] = DOMAIN_FLOOR_Y
            print(f"Particle missed the obstacle and hit the domain floor at y={DOMAIN_FLOOR_Y:.1f}.")
            is_falling = False # Stop simulation for this particle
            
        # time.sleep(0.05) # uncomment to slow down the printout

if __name__ == '__main__':
    # A particle that starts at x=0 will hit the obstacle
    simulate_fluid_particle(start_x=0, start_y=10)
    
    # A particle that starts at x=8 will miss the obstacle
    simulate_fluid_particle(start_x=8, start_y=10)
