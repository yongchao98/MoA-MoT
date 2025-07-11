import math

def solve_sonic_simulation():
    """
    Simulates the Sonic-like game to answer the two questions:
    1. How many coins does Sonic collect in total?
    2. At the apex of his jump, how many coins does Sonic collect on the same 
       vertical level as the highest coin he reaches?
    """
    # Constants and initial values from the problem description
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
    vel_x, vel_y = 0.0, 0.0

    # Generate coin positions
    # x-coordinates: 15, 17, ..., 45
    # y-coordinates: 11, 13, ..., 23
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    collected_coins = set()

    # The AABB intersection test checks if the distance between the centers
    # of two objects on each axis is less than half the sum of their sizes on that axis.
    # X-axis collision: abs(sonic_pos[0] - coin[0]) < (sonic_size[0] + coin_size[0]) / 2 => abs < 1.0
    # Y-axis collision: abs(sonic_pos[1] - coin[1]) < (sonic_size[1] + coin_size[1]) / 2 => abs < 1.5

    # Simulation loop
    for t in range(time_steps):
        # Update horizontal velocity
        vel_x += x_accel * delta_t
        
        # Apply jump impulse for a specific duration (t=16, 17, 18, 19)
        if 15 < t < 20:
            vel_y += jump_power * delta_t

        # Update position based on velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t

        # Apply gravity to vertical velocity
        vel_y += g * delta_t

        # Ground constraint: If Sonic falls below y=10, reset position and vertical velocity
        if sonic_pos[1] < 10:
            sonic_pos[1] = 10.0
            vel_y = 0.0

        # Check for collision with every coin in every frame
        for coin in coins:
            if abs(sonic_pos[0] - coin[0]) < 1.0 and abs(sonic_pos[1] - coin[1]) < 1.5:
                # Add coin to a set to ensure each coin is counted only once
                collected_coins.add(tuple(coin))

    # --- Part 1: Calculate the total number of unique coins collected ---
    total_collected_count = len(collected_coins)

    # --- Part 2: Calculate coins collected on the highest collected level ---
    highest_level_coin_count = 0
    if collected_coins:
        # First, find the highest y-coordinate among all the coins that were collected
        max_y_collected = max(coin[1] for coin in collected_coins)
        
        # Then, count how many of the collected coins share that y-coordinate
        highest_level_coin_count = sum(1 for coin in collected_coins if coin[1] == max_y_collected)
    
    # Print the final answer in the required format
    print(f"{total_collected_count} {highest_level_coin_count}")

# Run the simulation
solve_sonic_simulation()