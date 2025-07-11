# Here is a minimal Sonic-like game simulation.

def solve_sonic_coins():
    """
    Simulates a Sonic-like game to determine coin collection statistics.
    """
    # Game and simulation parameters
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5

    # Sonic's properties (using floating-point numbers for precision)
    sonic_pos, sonic_size = [5.0, 10.0], [1, 2]
    vel_x, vel_y = 0.0, 0.0

    # Generate the grid of coins
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]

    # --- Simulation Setup ---
    # Use a copy of the coins list to track remaining coins
    remaining_coins = list(coins)
    # Log to store the coordinates of all collected coins
    collected_coins_log = []

    def aabb_intersection(pos1, size1, pos2, size2):
        """Checks if two Axis-Aligned Bounding Boxes intersect."""
        x1, y1 = pos1
        w1, h1 = size1
        x2, y2 = pos2
        w2, h2 = size2
        
        # Check for overlap on x-axis and y-axis
        if (x1 < x2 + w2) and (x1 + w1 > x2) and \
           (y1 < y2 + h2) and (y1 + h1 > y2):
            return True
        return False

    # --- Main Simulation Loop ---
    for t in range(time_steps):
        # Update horizontal velocity and position
        vel_x += x_accel * delta_t
        sonic_pos[0] += vel_x * delta_t
        
        # Apply jump impulse for a short duration
        if 15 < t < 20:
            vel_y += jump_power * delta_t
            
        # Update vertical velocity (gravity) and position
        sonic_pos[1] += vel_y * delta_t
        vel_y += g * delta_t
        
        # Ground constraint
        if sonic_pos[1] < 10:
            sonic_pos[1] = 10.0
            vel_y = 0.0
            
        # --- Collision Detection ---
        newly_collected = []
        for coin_pos in remaining_coins:
            # AABB test: sonic vs coin. Coin size is [1, 1].
            if aabb_intersection(sonic_pos, sonic_size, coin_pos, [1, 1]):
                newly_collected.append(coin_pos)
                
        # Update logs if any coins were collected
        if newly_collected:
            collected_coins_log.extend(newly_collected)
            # This is a safe way to remove items while iterating
            remaining_coins = [c for c in remaining_coins if c not in newly_collected]

    # --- Post-Simulation Analysis ---

    # Part 1: How many coins does Sonic collect in total?
    total_coins_collected = len(collected_coins_log)

    # Part 2: How many coins are collected on the same vertical level as the highest coin reached?
    coins_at_highest_level = 0
    if total_coins_collected > 0:
        # Find the y-coordinate of the highest coin collected
        max_y = max(coin[1] for coin in collected_coins_log)
        
        # Count how many collected coins have this y-coordinate
        coins_at_highest_level = sum(1 for coin in collected_coins_log if coin[1] == max_y)

    # Print the final answers as two numbers separated by a space
    print(f"{total_coins_collected} {coins_at_highest_level}")

solve_sonic_coins()
<<<14 8>>>