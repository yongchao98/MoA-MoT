# Here is a minimal Sonic-like game simulation.

def solve_sonic_coins():
    """
    Simulates the Sonic-like game to determine coin collection statistics.
    """
    # Initial simulation parameters
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    
    # Use floating-point numbers for position and velocity for precision
    sonic_pos, sonic_size = [5.0, 10.0], [1, 2]
    vel_x, vel_y = 0.0, 0.0
    
    # Coin setup (each coin is size [1, 1])
    coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

    # Data structures for analysis
    collected_coins_set = set()
    collections_by_time = {} 
    max_sonic_y = -1.0
    apex_time = -1

    # Main simulation loop
    for t in range(time_steps):
        # Update Sonic's horizontal velocity and position
        vel_x += x_accel * delta_t
        sonic_pos[0] += vel_x * delta_t

        # Apply jump force during a specific time window
        if 15 < t < 20:
            vel_y += jump_power * delta_t
        
        # Update vertical position and apply gravity
        sonic_pos[1] += vel_y * delta_t
        vel_y += g * delta_t
        
        # Check for and resolve collision with the floor
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0
            
        # Track the apex (highest point) of the jump
        if sonic_pos[1] > max_sonic_y:
            max_sonic_y = sonic_pos[1]
            apex_time = t

        # Define Sonic's AABB for the current frame
        sonic_x1 = sonic_pos[0]
        sonic_y1 = sonic_pos[1]
        sonic_x2 = sonic_pos[0] + sonic_size[0]
        sonic_y2 = sonic_pos[1] + sonic_size[1]

        coins_collected_this_frame = []
        
        # Perform AABB intersection test against all coins
        for coin_pos in coins:
            coin_x1 = float(coin_pos[0])
            coin_y1 = float(coin_pos[1])
            coin_x2 = coin_x1 + 1.0
            coin_y2 = coin_y1 + 1.0

            # AABB overlap condition
            if sonic_x1 < coin_x2 and sonic_x2 > coin_x1 and \
               sonic_y1 < coin_y2 and sonic_y2 > coin_y1:
                
                # If a coin is collected, record it
                coin_tuple = tuple(coin_pos)
                collected_coins_set.add(coin_tuple)
                coins_collected_this_frame.append(coin_pos)
                
        if coins_collected_this_frame:
            collections_by_time[t] = coins_collected_this_frame

    # --- Post-simulation Analysis ---

    # 1. How many coins does Sonic collect in total?
    total_collected = len(collected_coins_set)

    # 2. At the apex of his jump, how many coins does Sonic collect 
    #    on the same vertical level as the highest coin he reaches?
    num_coins_at_highest_apex_level = 0
    coins_at_apex = collections_by_time.get(apex_time, [])

    if coins_at_apex:
        # Find the y-coordinate of the highest coin collected during the apex frame
        highest_coin_y_at_apex = max(coin[1] for coin in coins_at_apex)
        
        # Count how many coins from that frame share that highest y-coordinate
        num_coins_at_highest_apex_level = sum(1 for coin in coins_at_apex if coin[1] == highest_coin_y_at_apex)

    # Print the final answers separated by a space
    print(f"{total_collected} {num_coins_at_highest_apex_level}")

solve_sonic_coins()