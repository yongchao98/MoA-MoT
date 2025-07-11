def solve_sonic_coins():
    """
    Simulates a Sonic-like game to determine coin collection statistics.
    """
    # Here is a minimal Sonic-like game simulation.
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
    vel_x, vel_y = 0.0, 0.0
    initial_coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    remaining_coins = list(initial_coins)
    total_collected_coins = []
    collections_by_frame = {}

    max_sonic_y = -1.0
    apex_time_step = -1

    for t in range(time_steps):
        # Update velocity based on acceleration and jump
        vel_x += x_accel * delta_t
        if 15 < t < 20:
            vel_y += jump_power * delta_t

        # Update position based on velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t

        # Apply gravity to vertical velocity
        vel_y += g * delta_t

        # Ground collision check
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # Track the apex of the jump (highest point)
        if sonic_pos[1] > max_sonic_y:
            max_sonic_y = sonic_pos[1]
            apex_time_step = t

        # AABB intersection test for coin collection
        sonic_half_w, sonic_half_h = sonic_size[0] / 2.0, sonic_size[1] / 2.0
        sonic_x_min = sonic_pos[0] - sonic_half_w
        sonic_x_max = sonic_pos[0] + sonic_half_w
        sonic_y_min = sonic_pos[1] - sonic_half_h
        sonic_y_max = sonic_pos[1] + sonic_half_h

        coin_half_w, coin_half_h = coin_size[0] / 2.0, coin_size[1] / 2.0
        
        coins_to_keep = []
        collected_this_frame = []

        for coin_pos in remaining_coins:
            coin_x_min = coin_pos[0] - coin_half_w
            coin_x_max = coin_pos[0] + coin_half_w
            coin_y_min = coin_pos[1] - coin_half_h
            coin_y_max = coin_pos[1] + coin_half_h

            # Check for intersection
            if (sonic_x_min < coin_x_max and sonic_x_max > coin_x_min and
                sonic_y_min < coin_y_max and sonic_y_max > coin_y_min):
                collected_this_frame.append(coin_pos)
            else:
                coins_to_keep.append(coin_pos)
        
        if collected_this_frame:
            total_collected_coins.extend(collected_this_frame)
            collections_by_frame[t] = collected_this_frame
            remaining_coins = coins_to_keep

    # --- Analysis ---

    # 1. Total coins collected
    total_collected_count = len(total_collected_coins)

    # 2. Coins collected on the same vertical level as the highest coin at the apex
    coins_at_apex_frame = collections_by_frame.get(apex_time_step, [])
    
    apex_level_count = 0
    if coins_at_apex_frame:
        # Find the y-coordinate of the highest coin collected at the apex
        highest_coin_y_at_apex = -1.0
        for coin in coins_at_apex_frame:
            if coin[1] > highest_coin_y_at_apex:
                highest_coin_y_at_apex = coin[1]
        
        # Count how many coins at the apex share this highest y-level
        count = 0
        for coin in coins_at_apex_frame:
            if coin[1] == highest_coin_y_at_apex:
                count += 1
        apex_level_count = count

    # Print the final answers separated by a space
    print(f"{total_collected_count} {apex_level_count}")

solve_sonic_coins()