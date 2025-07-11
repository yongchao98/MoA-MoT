def solve_sonic_puzzle():
    """
    This function simulates the Sonic-like game to answer two questions:
    1. How many coins does Sonic collect in total?
    2. At the apex of his jump, how many coins does Sonic collect on the same vertical level as the highest coin he reaches?

    The simulation is run in two main stages. First, a full run to determine the total coins collected and the trajectory.
    Second, an analysis of the trajectory to find the apex of the jump, followed by a specific collision check at that moment.
    """
    
    # --- Parameters ---
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5
    
    # --- Initial State ---
    sonic_pos_init = [5.0, 10.0]
    sonic_size = [1.0, 2.0]
    
    # --- Coins ---
    all_coins = {(x, y) for x in range(15, 46, 2) for y in range(11, 24, 2)}
    
    # --- Full Simulation for Q1 and to record trajectory ---
    sonic_pos = list(sonic_pos_init)
    vel_x, vel_y = 0.0, 0.0
    collected_coins = set()
    trajectory_vel_y = {} # Store vertical velocity at the start of each time step
    
    for t in range(time_steps):
        # Store state at the beginning of the step (before updates)
        trajectory_vel_y[t] = vel_y
        
        # --- Physics Update ---
        vel_x += x_accel * delta_t
        if 15 < t < 20:  # Jump impulse applied at t = 16, 17, 18, 19
            vel_y += jump_power * delta_t
            
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t
        vel_y += g * delta_t
        
        # --- Floor Collision ---
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0
            
        # --- Coin Collection Check ---
        sonic_x_min = sonic_pos[0] - sonic_size[0] / 2.0
        sonic_x_max = sonic_pos[0] + sonic_size[0] / 2.0
        sonic_y_min = sonic_pos[1] - sonic_size[1] / 2.0
        sonic_y_max = sonic_pos[1] + sonic_size[1] / 2.0
        
        # Check against coins not yet collected to optimize
        coins_to_check = all_coins - collected_coins
        for coin in coins_to_check:
            cx, cy = coin
            coin_x_min, coin_x_max = cx - 0.5, cx + 0.5
            coin_y_min, coin_y_max = cy - 0.5, cy + 0.5
            
            # AABB intersection test
            if (sonic_x_max > coin_x_min and sonic_x_min < coin_x_max and
                sonic_y_max > coin_y_min and sonic_y_min < coin_y_max):
                collected_coins.add(coin)

    # --- Answer 1: Total coins collected ---
    answer1 = len(collected_coins)
    
    # --- Analysis for Q2 ---
    
    if not collected_coins:
        # If no coins were collected, the second answer is 0.
        print(f"{answer1} 0")
        return

    # Find the y-coordinate of the highest coin collected during the entire run.
    highest_collected_coin_y = max(c[1] for c in collected_coins)
    
    # Find the apex time step of the first jump.
    # The apex is the time step where the vertical velocity is closest to zero during the jump arc.
    apex_t = -1
    min_abs_vel_y = float('inf')
    
    # The jump starts at t=16. We need to find when it ends to define the search window.
    # The jump ends when Sonic returns to the ground (y=10) after being airborne.
    # We can determine this by re-simulating to find the landing time.
    _pos = list(sonic_pos_init)
    _vx, _vy = 0.0, 0.0
    landing_t = time_steps # Default to end of simulation if no landing
    for t_land in range(time_steps):
        _vx += x_accel * delta_t
        if 15 < t_land < 20: _vy += jump_power * delta_t
        _pos[0] += _vx * delta_t
        _pos[1] += _vy * delta_t
        _vy += g * delta_t
        if _pos[1] < 10.0:
            _pos[1], _vy = 10.0, 0.0
            if t_land > 20: # Ensure it's the landing after the jump impulse
                landing_t = t_land
                break
    
    # Find the apex within the jump arc (from t=16 to landing_t)
    for t in range(16, landing_t + 1):
        current_vel_y = trajectory_vel_y[t]
        if abs(current_vel_y) < min_abs_vel_y:
            min_abs_vel_y = abs(current_vel_y)
            apex_t = t
            
    # Re-simulate up to the apex time step to get the exact position for the collision check.
    sonic_pos_at_apex = list(sonic_pos_init)
    vel_x, vel_y = 0.0, 0.0
    
    for t in range(apex_t + 1):
        vel_x += x_accel * delta_t
        if 15 < t < 20:
            vel_y += jump_power * delta_t
            
        sonic_pos_at_apex[0] += vel_x * delta_t
        sonic_pos_at_apex[1] += vel_y * delta_t
        vel_y += g * delta_t
        
        if sonic_pos_at_apex[1] < 10.0:
            sonic_pos_at_apex[1], vel_y = 10.0, 0.0

    # Perform an AABB collision check at the apex position.
    sonic_x_min = sonic_pos_at_apex[0] - sonic_size[0] / 2.0
    sonic_x_max = sonic_pos_at_apex[0] + sonic_size[0] / 2.0
    sonic_y_min = sonic_pos_at_apex[1] - sonic_size[1] / 2.0
    sonic_y_max = sonic_pos_at_apex[1] + sonic_size[1] / 2.0
    
    coins_collected_at_apex = set()
    for coin in all_coins:
        cx, cy = coin
        coin_x_min, coin_x_max = cx - 0.5, cx + 0.5
        coin_y_min, coin_y_max = cy - 0.5, cy + 0.5
        
        if (sonic_x_max > coin_x_min and sonic_x_min < coin_x_max and
            sonic_y_max > coin_y_min and sonic_y_min < coin_y_max):
            coins_collected_at_apex.add(coin)
            
    # Count how many of the coins collected at the apex have the same y-level as the highest coin collected overall.
    answer2 = 0
    for coin in coins_collected_at_apex:
        if coin[1] == highest_collected_coin_y:
            answer2 += 1
            
    print(f"{answer1} {answer2}")

solve_sonic_puzzle()