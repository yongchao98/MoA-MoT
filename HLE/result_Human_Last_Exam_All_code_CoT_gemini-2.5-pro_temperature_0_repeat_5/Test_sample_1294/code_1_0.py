import math

# Here is a minimal Sonic-like game simulation.
# This script will simulate the game and calculate the required values.

def solve_sonic_coins():
    """
    This function encapsulates the entire simulation and calculation process.
    """

    # --- Helper function for AABB collision detection ---
    def check_collision(pos1, size1, pos2, size2):
        """
        Checks for intersection between two Axis-Aligned Bounding Boxes.
        Assumes pos is the bottom-left corner.
        """
        # Check for overlap on the X axis
        x_overlap = (pos1[0] < pos2[0] + size2[0] and 
                     pos1[0] + size1[0] > pos2[0])
        # Check for overlap on the Y axis
        y_overlap = (pos1[1] < pos2[1] + size2[1] and 
                     pos1[1] + size1[1] > pos2[1])
        return x_overlap and y_overlap

    # --- Part 1: Full Simulation to gather primary data ---

    # Initialize simulation parameters from the problem description
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    
    # Use floating-point numbers for positions and velocities for accuracy
    initial_sonic_pos = [5.0, 10.0]
    sonic_size = [1.0, 2.0]
    initial_vel_x, initial_vel_y = 0.0, 0.0
    
    # Generate coin positions
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    # State variables for the main simulation
    sonic_pos = list(initial_sonic_pos)
    vel_x, vel_y = initial_vel_x, initial_vel_y
    collected_coin_indices = set()
    max_sonic_y = -1.0
    apex_time = -1

    # Main simulation loop
    for t in range(time_steps):
        # Update horizontal velocity
        vel_x += x_accel * delta_t
        
        # Apply jump power during the specified time window
        if 15 < t < 20:
            vel_y += jump_power * delta_t
            
        # Update positions using current velocities
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        
        # Apply gravity to vertical velocity
        vel_y += g * delta_t
        
        # Check for and resolve ground collision
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # Track the highest point (apex) of Sonic's jump
        if sonic_pos[1] > max_sonic_y:
            max_sonic_y = sonic_pos[1]
            apex_time = t

        # Check for collisions with all coins
        for i, coin_pos in enumerate(coins):
            if i not in collected_coin_indices:
                if check_collision(sonic_pos, sonic_size, coin_pos, coin_size):
                    collected_coin_indices.add(i)

    # --- Analysis and Part 2 Calculation ---

    # Answer 1: Total number of unique coins collected
    total_collected = len(collected_coin_indices)

    # Determine the y-coordinate of the highest coin Sonic collected
    max_collected_coin_y = -1.0
    if collected_coin_indices:
        # Find the coin with the maximum y-coordinate from the collected set
        highest_collected_coin = max((coins[i] for i in collected_coin_indices), key=lambda c: c[1])
        max_collected_coin_y = highest_collected_coin[1]

    # To find the coins collected *at* the apex, we need Sonic's exact position at that time.
    # We re-simulate up to the apex_time to get this state.
    
    # Reset Sonic to initial state
    sonic_pos = list(initial_sonic_pos)
    vel_x, vel_y = initial_vel_x, initial_vel_y
    
    apex_sonic_pos = list(sonic_pos)
    if apex_time > -1:
        for t in range(apex_time + 1):
            vel_x += x_accel * delta_t
            if 15 < t < 20:
                vel_y += jump_power * delta_t
            sonic_pos[0] += vel_x * delta_t
            sonic_pos[1] += vel_y * delta_t
            vel_y += g * delta_t
            if sonic_pos[1] < 10.0:
                sonic_pos[1], vel_y = 10.0, 0.0
        apex_sonic_pos = list(sonic_pos)

    # Answer 2: At the apex, count how many coins are collected on the highest level
    apex_coins_on_level = 0
    if apex_time > -1 and max_collected_coin_y > -1.0:
        for coin_pos in coins:
            # Check for collision at the apex position
            if check_collision(apex_sonic_pos, sonic_size, coin_pos, coin_size):
                # If it collides, check if it's on the target y-level
                if math.isclose(coin_pos[1], max_collected_coin_y):
                    apex_coins_on_level += 1

    # Print the final answers as two numbers separated by a space
    print(f"{total_collected} {apex_coins_on_level}")

# Execute the function to print the results
solve_sonic_coins()