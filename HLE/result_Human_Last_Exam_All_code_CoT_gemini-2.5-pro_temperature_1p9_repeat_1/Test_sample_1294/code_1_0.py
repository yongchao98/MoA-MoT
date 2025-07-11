import math

# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
# Use floating point numbers for position for better precision
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
# Coins are also defined with float coordinates
coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# A set to store the coordinates of all unique coins collected
collected_coins_set = set()
# A dictionary to store which coins were collected at each time step
collections_per_frame = {}
# Variables to track the apex (highest point) of the jump
max_sonic_y = -1.0
apex_time = -1

# Main simulation loop
for t in range(time_steps):
    # Update horizontal velocity
    vel_x += x_accel * delta_t
    # Apply jump power during the specified window
    if 15 < t < 20:
        vel_y += jump_power * delta_t
        
    # Update Sonic's position based on velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to vertical velocity
    vel_y += g * delta_t
    
    # Ground check: if Sonic is below the ground, reset his position and vertical velocity
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0

    # Track the apex by finding the maximum y-position and the time it occurred
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]
        apex_time = t

    # AABB collision check logic
    sonic_half_w = sonic_size[0] / 2.0
    sonic_half_h = sonic_size[1] / 2.0
    # Coin size is [1, 1], so its half-width/height is 0.5
    coin_half_size = 0.5
    
    # A temporary set to hold coins collected in this frame
    newly_collected_this_frame = set()

    # Iterate through all coins to check for collisions
    for coin_pos in coins:
        coin_tuple = tuple(coin_pos)
        # Skip coins that have already been collected
        if coin_tuple in collected_coins_set:
            continue
            
        # AABB intersection test: check if collision boxes overlap on both X and Y axes
        is_x_overlapping = abs(sonic_pos[0] - coin_pos[0]) < (sonic_half_w + coin_half_size)
        is_y_overlapping = abs(sonic_pos[1] - coin_pos[1]) < (sonic_half_h + coin_half_size)
        
        if is_x_overlapping and is_y_overlapping:
            collected_coins_set.add(coin_tuple)
            newly_collected_this_frame.add(coin_tuple)

    # If any coins were collected, log them for this time step
    if newly_collected_this_frame:
        collections_per_frame[t] = newly_collected_this_frame

# --- Question 1: How many coins does Sonic collect in total? ---
total_collected_count = len(collected_coins_set)

# --- Question 2: At the apex of his jump, how many coins does Sonic collect 
# --- on the same vertical level as the highest coin he reaches? ---
coins_at_apex_level = 0
if collected_coins_set:
    # First, find the y-level of the highest coin collected overall
    highest_coin_y_overall = max(c[1] for c in collected_coins_set)
    
    # Get the set of coins collected during the apex frame (if any)
    apex_frame_collections = collections_per_frame.get(apex_time, set())
    
    if apex_frame_collections:
        # Count how many coins from the apex frame collection are at the highest level
        coins_at_apex_level = sum(1 for coin in apex_frame_collections if coin[1] == highest_coin_y_overall)

# Print the final two numbers separated by a space
print(f"{total_collected_count} {coins_at_apex_level}")