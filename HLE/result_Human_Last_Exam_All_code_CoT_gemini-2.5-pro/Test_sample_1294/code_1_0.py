# Here is a minimal Sonic-like game simulation.

def solve_sonic_coins():
    """
    Simulates Sonic's movement to determine coin collection statistics.
    """
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    # Use floating point numbers for positions to maintain precision
    sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
    vel_x, vel_y = 0.0, 0.0
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    # Use a set to store indices of collected coins to prevent duplicates
    collected_coin_indices = set()

    for t in range(time_steps):
        # Update velocities
        vel_x += x_accel * delta_t
        if 15 < t < 20:  # Jump is applied during steps 16, 17, 18, 19
            vel_y += jump_power * delta_t
        
        # Apply gravity
        vel_y += g * delta_t
        
        # Update position
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        
        # Ground constraint
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # Define Sonic's AABB for collision detection
        sonic_left = sonic_pos[0]
        sonic_right = sonic_pos[0] + sonic_size[0]
        sonic_bottom = sonic_pos[1]
        sonic_top = sonic_pos[1] + sonic_size[1]

        # Check for collision with each coin
        for i, coin_pos in enumerate(coins):
            # Skip coins that have already been collected
            if i in collected_coin_indices:
                continue

            # Define coin's AABB
            coin_left = coin_pos[0]
            coin_right = coin_pos[0] + coin_size[0]
            coin_bottom = coin_pos[1]
            coin_top = coin_pos[1] + coin_size[1]
            
            # AABB intersection test
            is_colliding = (sonic_right > coin_left and 
                            sonic_left < coin_right and 
                            sonic_top > coin_bottom and 
                            sonic_bottom < coin_top)

            if is_colliding:
                collected_coin_indices.add(i)

    # --- Part 1: Calculate the total number of collected coins ---
    total_collected = len(collected_coin_indices)

    # --- Part 2: Calculate coins collected at the highest vertical level ---
    count_at_highest_level = 0
    if total_collected > 0:
        # Get the positions of the coins that were actually collected
        collected_coins_list = [coins[i] for i in collected_coin_indices]
        
        # Find the maximum y-coordinate among the collected coins
        highest_coin_y = max(coin[1] for coin in collected_coins_list)
        
        # Count how many collected coins are at that highest y-level
        count_at_highest_level = sum(1 for coin in collected_coins_list if coin[1] == highest_coin_y)

    # Print the two answers separated by a space
    print(f"{total_collected} {count_at_highest_level}")

solve_sonic_coins()
<<<23 2>>>