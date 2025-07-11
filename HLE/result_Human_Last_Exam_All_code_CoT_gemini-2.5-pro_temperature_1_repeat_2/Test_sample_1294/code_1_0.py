import math

# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# List to keep track of all coins that have been collected
collected_coins_log = []

# Main simulation loop
for t in range(time_steps):
    # Update horizontal velocity and apply jump power
    vel_x += x_accel * delta_t
    if 15 < t < 20:
        vel_y += jump_power * delta_t
        
    # Update Sonic's position based on velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to vertical velocity
    vel_y += g * delta_t
    
    # Check for and handle collision with the ground (at y=10)
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0

    # --- Coin Collection Logic ---
    # Define Sonic's bounding box for the current frame
    sonic_x_min = sonic_pos[0]
    sonic_x_max = sonic_pos[0] + sonic_size[0]
    sonic_y_min = sonic_pos[1]
    sonic_y_max = sonic_pos[1] + sonic_size[1]
    
    # A temporary list to hold coins collected in this frame
    coins_collected_this_frame = []
    
    # Check for collision with each coin
    for coin in coins:
        # Define the coin's bounding box (size is [1, 1])
        coin_x_min = coin[0]
        coin_x_max = coin[0] + 1.0
        coin_y_min = coin[1]
        coin_y_max = coin[1] + 1.0
        
        # AABB intersection test
        is_overlapping_x = (sonic_x_min < coin_x_max and sonic_x_max > coin_x_min)
        is_overlapping_y = (sonic_y_min < coin_y_max and sonic_y_max > coin_y_min)
        
        if is_overlapping_x and is_overlapping_y:
            coins_collected_this_frame.append(coin)
            
    # Process the collected coins for this frame
    if coins_collected_this_frame:
        for coin in coins_collected_this_frame:
            collected_coins_log.append(coin)
            coins.remove(coin) # Remove the coin so it can't be collected again

# --- Analysis of Results ---

# Part 1: How many coins does Sonic collect in total?
total_coins_collected = len(collected_coins_log)

# Part 2: At the apex of his jump, how many coins does Sonic collect 
# on the same vertical level as the highest coin he reaches?
apex_level_coins = 0
if total_coins_collected > 0:
    # Find the maximum y-coordinate among all collected coins
    highest_coin_y_level = max(coin[1] for coin in collected_coins_log)
    
    # Count how many of the collected coins are on that same highest level
    apex_level_coins = sum(1 for coin in collected_coins_log if coin[1] == highest_coin_y_level)

# Print the final answers separated by a space
print(total_coins_collected, apex_level_coins)
<<<21 7>>>