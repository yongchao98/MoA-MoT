import collections

# Here is a minimal Sonic-like game simulation.
# The original code is augmented to track collections and find the answers.

# --- Initial Simulation Setup ---
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
# Using floats for position for better accuracy
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
# Coin size is [1, 1]
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# --- Tracking Variables ---
# Use a set to store unique collected coins
collected_coins_total = set()
# Store coins collected at each specific time step
collections_by_frame = collections.defaultdict(set)
# Track the apex of the jump
max_sonic_y = sonic_pos[1]
apex_time_step = 0

# --- Main Simulation Loop ---
for t in range(time_steps):
    # Update horizontal velocity and position
    vel_x += x_accel * delta_t
    sonic_pos[0] += vel_x * delta_t
    
    # Apply jump power during the specified interval
    if 15 < t < 20: # For t = 16, 17, 18, 19
        vel_y += jump_power * delta_t
    
    # Update vertical velocity and position with gravity
    vel_y += g * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Check for and resolve collision with the ground
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0
        
    # Check if a new maximum height (apex) has been reached
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]
        apex_time_step = t
        
    # --- AABB Intersection Test for Coin Collection ---
    sonic_x_min = sonic_pos[0] - sonic_size[0] / 2
    sonic_x_max = sonic_pos[0] + sonic_size[0] / 2
    sonic_y_min = sonic_pos[1] - sonic_size[1] / 2
    sonic_y_max = sonic_pos[1] + sonic_size[1] / 2
    
    coin_half_size = 0.5
    for coin in coins:
        coin_x, coin_y = coin[0], coin[1]
        coin_x_min = coin_x - coin_half_size
        coin_x_max = coin_x + coin_half_size
        coin_y_min = coin_y - coin_half_size
        coin_y_max = coin_y + coin_half_size
        
        # Check for overlap on both axes
        if (sonic_x_max > coin_x_min and sonic_x_min < coin_x_max and
            sonic_y_max > coin_y_min and sonic_y_min < coin_y_max):
            
            coin_tuple = (coin_x, coin_y)
            collected_coins_total.add(coin_tuple)
            collections_by_frame[t].add(coin_tuple)

# --- Analysis and Final Answer Calculation ---

# 1. How many coins does Sonic collect in total?
total_collected_count = len(collected_coins_total)

# 2. At the apex of his jump, how many coins does Sonic collect on the same vertical level as the highest coin he reaches?

# Find the highest y-coordinate among all collected coins
if not collected_coins_total:
    highest_collected_coin_y = -1
else:
    highest_collected_coin_y = max(c[1] for c in collected_coins_total)

# Get the coins collected during the apex frame
coins_at_apex = collections_by_frame[apex_time_step]

# Count how many of those coins are on the highest level
count_at_apex_on_highest_level = 0
if highest_collected_coin_y != -1:
    for coin in coins_at_apex:
        if coin[1] == highest_collected_coin_y:
            count_at_apex_on_highest_level += 1

# Print the final answer
print(f"{total_collected_count} {count_at_apex_on_highest_level}")