def solve():
    """
    Here is a minimal Sonic-like game simulation.
    """
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    sonic_pos, sonic_size = [5.0, 10.0], [1, 2]  # Use floats for position
    vel_x, vel_y = 0.0, 0.0
    coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1, 1]

    # Use a set to store the coordinates of collected coins to avoid duplicates
    collected_coins = set()

    for t in range(time_steps):
        # Update horizontal velocity
        vel_x += x_accel * delta_t
        
        # Apply jump power if within the jump time frame
        if 15 < t < 20:
            vel_y += jump_power * delta_t
            
        # Update position based on velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        
        # Apply gravity to vertical velocity
        vel_y += g * delta_t
        
        # Check for ground collision
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # AABB Collision Detection
        sonic_x_min = sonic_pos[0]
        sonic_x_max = sonic_pos[0] + sonic_size[0]
        sonic_y_min = sonic_pos[1]
        sonic_y_max = sonic_pos[1] + sonic_size[1]

        for coin in coins:
            coin_x_min = coin[0]
            coin_x_max = coin[0] + coin_size[0]
            coin_y_min = coin[1]
            coin_y_max = coin[1] + coin_size[1]

            # Check for overlap on both axes
            if (sonic_x_min < coin_x_max and sonic_x_max > coin_x_min and
                sonic_y_min < coin_y_max and sonic_y_max > coin_y_min):
                collected_coins.add(tuple(coin))

    # --- Part 1: Calculate total coins collected ---
    total_collected = len(collected_coins)

    # --- Part 2: Calculate coins collected at the highest level reached ---
    apex_level_count = 0
    if collected_coins:
        # Find the maximum y-coordinate among all collected coins
        y_highest_coin = -1
        for coin in collected_coins:
            if coin[1] > y_highest_coin:
                y_highest_coin = coin[1]
        
        # Count how many of the collected coins are on that same vertical level
        count = 0
        for coin in collected_coins:
            if coin[1] == y_highest_coin:
                count += 1
        apex_level_count = count

    # Print the final answer
    print(f"{total_collected} {apex_level_count}")

solve()
<<<18 2>>>