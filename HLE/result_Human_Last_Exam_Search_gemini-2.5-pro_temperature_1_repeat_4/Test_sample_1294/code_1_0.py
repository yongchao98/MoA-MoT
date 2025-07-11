def solve_sonic_coins():
    """
    Simulates the Sonic-like game to answer two questions:
    1. How many coins does Sonic collect in total?
    2. At the apex of his jump, how many coins does Sonic collect on the same vertical level as the highest coin he reaches?
    """

    # Simulation parameters
    delta_t = 0.05
    time_steps = 1000
    g = -9.8
    jump_power = 75
    x_accel = 5.5

    # Sonic's properties
    sonic_pos = [5.0, 10.0]
    sonic_size = [1.0, 2.0]
    vel_x = 0.0
    vel_y = 0.0

    # Coin properties
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    collected_coins_indices = set()
    max_y = -1.0
    apex_info = {'t': -1, 'pos': [0,0], 'collected_coins_at_apex': []}

    # Simulation loop
    for t in range(time_steps):
        # Update horizontal velocity
        vel_x += x_accel * delta_t

        # Apply jump force
        if 15 < t < 20:
            vel_y += jump_power * delta_t

        # Update position
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t

        # Apply gravity
        vel_y += g * delta_t

        # Floor collision
        if sonic_pos[1] < 10.0:
            sonic_pos[1] = 10.0
            vel_y = 0.0
        
        # Check for apex
        if sonic_pos[1] > max_y:
            max_y = sonic_pos[1]
            apex_info['t'] = t
            apex_info['pos'] = list(sonic_pos)
            apex_info['collected_coins_at_apex'] = []


        # AABB collision detection
        sonic_min_x = sonic_pos[0]
        sonic_max_x = sonic_pos[0] + sonic_size[0]
        sonic_min_y = sonic_pos[1]
        sonic_max_y = sonic_pos[1] + sonic_size[1]
        
        coins_collected_this_frame = []

        for i, coin_pos in enumerate(coins):
            coin_min_x = coin_pos[0]
            coin_max_x = coin_pos[0] + coin_size[0]
            coin_min_y = coin_pos[1]
            coin_max_y = coin_pos[1] + coin_size[1]

            if (sonic_min_x < coin_max_x and sonic_max_x > coin_min_x and
                    sonic_min_y < coin_max_y and sonic_max_y > coin_min_y):
                collected_coins_indices.add(i)
                coins_collected_this_frame.append(coin_pos)

        if t == apex_info['t']:
             apex_info['collected_coins_at_apex'].extend(coins_collected_this_frame)


    # Part 1: Total collected coins
    total_collected_coins = len(collected_coins_indices)

    # Part 2: Coins collected at the apex
    highest_coin_y_at_apex = 0
    if apex_info['collected_coins_at_apex']:
        for coin in apex_info['collected_coins_at_apex']:
            if coin[1] > highest_coin_y_at_apex:
                highest_coin_y_at_apex = coin[1]

    
    count_at_apex_level = 0
    if apex_info['collected_coins_at_apex']:
        for coin in apex_info['collected_coins_at_apex']:
            if coin[1] == highest_coin_y_at_apex:
                count_at_apex_level += 1


    return total_collected_coins, count_at_apex_level

total_coins, apex_coins = solve_sonic_coins()
print(f"{total_coins} {apex_coins}")