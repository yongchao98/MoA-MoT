import math

def solve():
    """
    Simulates the Sonic-like game to determine coin collection statistics.
    """
    # Here is a minimal Sonic-like game simulation.
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    # Use floating-point numbers for positions and velocities for precision
    sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
    vel_x, vel_y = 0.0, 0.0
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    
    # --- Simulation setup ---
    coin_size = [1.0, 1.0]
    collected_coins_log = []

    # Pre-calculate collision distance checks based on AABB test
    # abs(pos1 - pos2) < (size1/2 + size2/2)
    x_dist_check = (sonic_size[0] + coin_size[0]) / 2.0
    y_dist_check = (sonic_size[1] + coin_size[1]) / 2.0

    # --- Main simulation loop ---
    for t in range(time_steps):
        # Update horizontal velocity and apply jump impulse
        vel_x += x_accel * delta_t
        if 15 < t < 20: # True for t = 16, 17, 18, 19
            vel_y += jump_power * delta_t
        
        # Update Sonic's position based on velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        
        # Apply gravity to vertical velocity
        vel_y += g * delta_t
        
        # Check for ground collision and reset state if necessary
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # --- Collision detection and collection ---
        remaining_coins = []
        collected_this_frame = False
        for coin in coins:
            # Check for AABB collision
            if abs(sonic_pos[0] - coin[0]) < x_dist_check and \
               abs(sonic_pos[1] - coin[1]) < y_dist_check:
                
                # Log the collected coin's data
                collected_coins_log.append({'pos': coin, 't': t})
                collected_this_frame = True
            else:
                # If not collected, the coin remains for the next frame
                remaining_coins.append(coin)
        
        # Update the list of active coins
        if collected_this_frame:
            coins = remaining_coins

    # --- Analysis of collected data ---

    # Question 1: How many coins does Sonic collect in total?
    total_collected_coins = len(collected_coins_log)

    # Question 2: Coins collected at the apex.
    apex_level_coins_count = 0
    if total_collected_coins > 0:
        # Find the highest y-coordinate of any collected coin
        max_y_collected = max(info['pos'][1] for info in collected_coins_log)

        # Find the earliest timestep when a coin at this max_y was collected.
        # This defines the "apex" frame for the purpose of the question.
        apex_time = min(info['t'] for info in collected_coins_log if info['pos'][1] == max_y_collected)

        # Count how many coins were collected at that specific time (apex_time)
        # AND are on the same vertical level (max_y_collected).
        count = 0
        for info in collected_coins_log:
            if info['t'] == apex_time and info['pos'][1] == max_y_collected:
                count += 1
        apex_level_coins_count = count
    
    # Print the final answer as two numbers separated by a space.
    print(f"{total_collected_coins} {apex_level_coins_count}")

solve()
<<<91 1>>>