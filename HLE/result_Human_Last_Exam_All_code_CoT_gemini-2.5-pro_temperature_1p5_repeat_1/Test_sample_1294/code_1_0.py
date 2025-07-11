# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5.0, 10.0], [1, 2]
vel_x, vel_y = 0.0, 0.0
coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# Use a set to store collected coins to avoid duplicates
collected_coins_set = set()

# Variables to find Sonic's state at the apex of his jump
max_sonic_y = -1.0
sonic_pos_at_apex = [0.0, 0.0]

# --- Main simulation loop ---
for t in range(time_steps):
    # Update velocity and position based on physics
    vel_x += x_accel * delta_t
    
    # Apply jump power during a specific time window
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    vel_y += g * delta_t
    
    # Ground check
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0

    # Check for the apex of the jump (highest point reached)
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]
        sonic_pos_at_apex = list(sonic_pos) # Store a copy of the position

    # --- AABB collision detection for collecting coins ---
    s_x, s_y = sonic_pos[0], sonic_pos[1]
    s_w, s_h = sonic_size[0], sonic_size[1]
    
    for coin in coins:
        # Skip coins that have already been collected
        if tuple(coin) in collected_coins_set:
            continue
            
        c_x, c_y = coin[0], coin[1]
        c_w, c_h = 1.0, 1.0  # Coin size is [1,1]

        # AABB intersection test conditions:
        # Overlap in x: sonic_x < coin_x + coin_width AND sonic_x + sonic_width > coin_x
        # Overlap in y: sonic_y < coin_y + coin_height AND sonic_y + sonic_height > coin_y
        if (s_x < c_x + c_w and s_x + s_w > c_x and
            s_y < c_y + c_h and s_y + s_h > c_y):
            collected_coins_set.add(tuple(coin))

# --- Post-simulation analysis and results ---

# Question 1: How many coins does Sonic collect in total?
total_collected_coins = len(collected_coins_set)

# Question 2: At the apex of his jump, how many coins does Sonic collect
# on the same vertical level as the highest coin he reaches?
apex_coins_count = 0

if total_collected_coins > 0:
    # Find the y-level of the highest coin that was actually collected
    y_highest_coin = max(c[1] for c in collected_coins_set)
    
    # Get all coins from the original grid that exist on this highest collected level
    coins_on_highest_level = [c for c in coins if c[1] == y_highest_coin]
    
    # Get Sonic's hitbox at the exact moment of the apex
    s_x_apex, s_y_apex = sonic_pos_at_apex[0], sonic_pos_at_apex[1]
    s_w_apex, s_h_apex = sonic_size[0], sonic_size[1]
    
    # Check for collision with these specific coins using Sonic's apex position
    for coin in coins_on_highest_level:
        c_x, c_y = coin[0], coin[1]
        c_w, c_h = 1.0, 1.0
        
        # AABB test at the apex
        if (s_x_apex < c_x + c_w and s_x_apex + s_w_apex > c_x and
            s_y_apex < c_y + c_h and s_y_apex + s_h_apex > c_y):
            apex_coins_count += 1

# Print the final answers
print(f"{total_collected_coins} {apex_coins_count}")
<<<45 2>>>