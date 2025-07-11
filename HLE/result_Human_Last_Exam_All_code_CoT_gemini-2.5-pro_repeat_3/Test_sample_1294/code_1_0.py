
# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
# Use floating point numbers for positions for better precision
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# A set to store the coordinates of collected coins to prevent duplicates
collected_coins = set()

# AABB collision parameters
sonic_width, sonic_height = sonic_size[0], sonic_size[1]
coin_width, coin_height = 1.0, 1.0
# Pre-calculate the distance thresholds for collision
x_dist_threshold = (sonic_width + coin_width) / 2
y_dist_threshold = (sonic_height + coin_height) / 2

# Main simulation loop
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
    
    # Check for ground collision
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0

    # Check for coin collection
    for coin in coins:
        coin_tuple = tuple(coin)
        # Skip coins that have already been collected
        if coin_tuple in collected_coins:
            continue
        
        # AABB intersection test: check if the distance between centers is
        # less than half the sum of the sizes on each axis.
        is_x_overlap = abs(sonic_pos[0] - coin[0]) < x_dist_threshold
        is_y_overlap = abs(sonic_pos[1] - coin[1]) < y_dist_threshold
        
        if is_x_overlap and is_y_overlap:
            collected_coins.add(coin_tuple)

# --- Post-simulation analysis ---

# 1. How many coins does Sonic collect in total?
total_collected = len(collected_coins)

# 2. How many coins does Sonic collect on the same vertical level
#    as the highest coin he reaches?

# Find the y-coordinate of the highest coin that was collected.
highest_coin_y = 0
if collected_coins:
    highest_coin_y = max(coin[1] for coin in collected_coins)

# Count how many of the collected coins are at that highest y-level.
apex_level_count = 0
if highest_coin_y > 0:
    apex_level_count = sum(1 for coin in collected_coins if coin[1] == highest_coin_y)

# Print the final answers as two numbers separated by a space.
print(total_collected, apex_level_count)
