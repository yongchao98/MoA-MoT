import math

def solve_sonic_coins():
    """
    Simulates a Sonic-like game to determine coin collection statistics.
    1. Calculates the total number of coins Sonic collects.
    2. At the apex of Sonic's jump, identifies the highest vertical level of a collected coin
       and counts how many coins were collected on that same level in that single frame.
    """
    # Initialize game parameters and state variables
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    # Use float lists for precise position and velocity calculations
    sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
    vel_x, vel_y = 0.0, 0.0
    # Generate coin positions as floats
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    # Variables for tracking results
    total_collected_coins = 0
    # history will store {'t': t, 'pos': [x,y], 'collected': [list_of_coins]} for each frame
    history = []

    # Helper function for AABB (Axis-Aligned Bounding Box) collision detection
    def check_collision(pos1, size1, pos2, size2):
        x1, y1 = pos1
        w1, h1 = size1
        x2, y2 = pos2
        w2, h2 = size2

        # Calculate bounding box edges for both objects
        x1_min, x1_max = x1 - w1 / 2, x1 + w1 / 2
        y1_min, y1_max = y1 - h1 / 2, y1 + h1 / 2
        x2_min, x2_max = x2 - w2 / 2, x2 + w2 / 2
        y2_min, y2_max = y2 - h2 / 2, y2 + h2 / 2

        # Return True if boxes overlap, False otherwise
        return x1_min < x2_max and x1_max > x2_min and y1_min < y2_max and y1_max > y2_min

    # Main simulation loop
    for t in range(time_steps):
        # Update horizontal velocity and apply jump force
        vel_x += x_accel * delta_t
        if 15 < t < 20:
            vel_y += jump_power * delta_t
        
        # Update position based on velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        
        # Apply gravity
        vel_y += g * delta_t
        
        # Apply ground constraint
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # Check for coin collisions
        collected_this_frame = []
        remaining_coins = []
        for coin in coins:
            if check_collision(sonic_pos, sonic_size, coin, coin_size):
                collected_this_frame.append(coin)
            else:
                remaining_coins.append(coin)
        
        # Update total count and the list of available coins
        if collected_this_frame:
            total_collected_coins += len(collected_this_frame)
            coins = remaining_coins

        # Record the state of this frame for later analysis
        history.append({
            't': t,
            'pos': list(sonic_pos),  # Store a copy of the position
            'collected': collected_this_frame
        })

    # --- Post-Simulation Analysis ---

    # Part 1: Total coins collected is already calculated
    total_coins_answer = total_collected_coins

    # Part 2: Find coins collected at the apex of the jump
    # Find the frame where Sonic's y-position was maximum
    apex_frame = max(history, key=lambda frame: frame['pos'][1])
    
    # Get the list of coins collected in that specific frame
    collections_at_apex = apex_frame['collected']
    
    apex_coins_answer = 0
    if collections_at_apex:
        # Find the maximum y-coordinate among the coins collected at the apex
        max_y_at_apex = max(coin[1] for coin in collections_at_apex)
        
        # Count how many of those coins were on that same vertical level
        apex_coins_answer = sum(1 for coin in collections_at_apex if coin[1] == max_y_at_apex)

    # Print the final answers
    print(f"{total_coins_answer} {apex_coins_answer}")

solve_sonic_coins()
<<<88 3>>>