# Here is a minimal Sonic-like game simulation.
def solve_sonic_coins():
    """
    Simulates a Sonic-like game to determine coin collection statistics.
    """
    # Initialize simulation parameters and game world state
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    sonic_pos, sonic_size = [5, 10], [1, 2]
    vel_x, vel_y = 0, 0
    
    # Generate coins and store them in a set of tuples for efficient operations
    coins_list = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coins_set = {tuple(c) for c in coins_list}
    
    # Use a set to store collected coins, automatically handling duplicates
    collected_coins_set = set()

    # Main simulation loop
    for t in range(time_steps):
        # Update horizontal velocity
        vel_x += x_accel * delta_t
        
        # Apply jump power during the specified time window
        if 15 < t < 20: # Corresponds to t = 16, 17, 18, 19
            vel_y += jump_power * delta_t
            
        # Update Sonic's position based on velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        
        # Apply gravity to vertical velocity
        vel_y += g * delta_t
        
        # Ground constraint: if Sonic is below the ground, reset his position and vertical velocity
        if sonic_pos[1] < 10:
            sonic_pos[1], vel_y = 10, 0
            
        # Define Sonic's AABB (Axis-Aligned Bounding Box) for collision detection
        sonic_x1 = sonic_pos[0]
        sonic_y1 = sonic_pos[1]
        sonic_x2 = sonic_pos[0] + sonic_size[0]
        sonic_y2 = sonic_pos[1] + sonic_size[1]
        
        # Identify coins collected in the current frame
        coins_to_collect_this_frame = set()
        for coin_pos_tuple in coins_set:
            coin_x1 = coin_pos_tuple[0]
            coin_y1 = coin_pos_tuple[1]
            # Coins have a size of [1, 1]
            coin_x2 = coin_pos_tuple[0] + 1
            coin_y2 = coin_pos_tuple[1] + 1
            
            # AABB intersection test to check for overlap
            if (sonic_x1 < coin_x2 and sonic_x2 > coin_x1 and
                sonic_y1 < coin_y2 and sonic_y2 > coin_y1):
                coins_to_collect_this_frame.add(coin_pos_tuple)
        
        # Add collected coins to the main set and remove them from the game world
        if coins_to_collect_this_frame:
            collected_coins_set.update(coins_to_collect_this_frame)
            coins_set.difference_update(coins_to_collect_this_frame)

    # --- Analysis of collected coins after simulation ---

    # Question 1: How many coins does Sonic collect in total?
    total_coins_collected = len(collected_coins_set)
    
    # Question 2: At the apex of his jump, how many coins does Sonic collect
    # on the same vertical level as the highest coin he reaches?
    # This is interpreted as: of all collected coins, what is the count of
    # coins at the highest y-level achieved?
    
    highest_level_coin_count = 0
    if collected_coins_set:
        # Find the maximum y-coordinate among all collected coins
        max_y = max(c[1] for c in collected_coins_set)
        
        # Count how many of the collected coins are at this maximum y-level
        highest_level_coin_count = sum(1 for c in collected_coins_set if c[1] == max_y)

    # Print the final answers
    print(f"{total_coins_collected} {highest_level_coin_count}")

solve_sonic_coins()
<<<36 4>>>