import collections

# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0

# Use a set of tuples for coins for efficient removal and lookup.
coins = set((float(x), float(y)) for x in range(15, 46, 2) for y in range(11, 24, 2))

# --- Data structures for tracking results ---
# Use a set to store unique collected coins
collected_coins_total = set()
# Variables to track the apex of the jump
max_sonic_y = -1.0
coins_collected_at_apex_frame = []

def check_aabb_collision(pos1, size1, pos2, size2):
    """Checks for AABB collision between two objects."""
    half_w1, half_h1 = size1[0] / 2, size1[1] / 2
    half_w2, half_h2 = size2[0] / 2, size2[1] / 2

    # Object 1's boundaries
    left1, right1 = pos1[0] - half_w1, pos1[0] + half_w1
    bottom1, top1 = pos1[1] - half_h1, pos1[1] + half_h1

    # Object 2's boundaries
    left2, right2 = pos2[0] - half_w2, pos2[0] + half_w2
    bottom2, top2 = pos2[1] - half_h2, pos2[1] + half_h2

    # Intersection test
    return left1 < right2 and right1 > left2 and bottom1 < top2 and top1 > bottom2

# --- Main simulation loop ---
for t in range(time_steps):
    # Update velocity based on acceleration and jump input
    vel_x += x_accel * delta_t
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    # Update position based on velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity
    vel_y += g * delta_t

    # Ground constraint
    if sonic_pos[1] < 10:
        sonic_pos[1], vel_y = 10.0, 0.0

    # --- Apex tracking logic ---
    # If we've reached a new highest point, this is our new potential apex frame
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]
        coins_collected_at_apex_frame = [] # Reset for the new apex

    # --- Collision detection and collection ---
    coins_to_remove_this_frame = set()
    coins_collected_this_frame = []

    for coin_pos in coins:
        if check_aabb_collision(sonic_pos, sonic_size, list(coin_pos), [1.0, 1.0]):
            collected_coins_total.add(coin_pos)
            coins_to_remove_this_frame.add(coin_pos)
            coins_collected_this_frame.append(coin_pos)
    
    # Remove the collected coins from the main set so they can't be collected again
    if coins_to_remove_this_frame:
        coins.difference_update(coins_to_remove_this_frame)

    # If this frame is the current apex, store the coins collected
    if sonic_pos[1] == max_sonic_y:
        coins_collected_at_apex_frame.extend(coins_collected_this_frame)

# --- Post-simulation analysis ---

# 1. How many coins does Sonic collect in total?
total_collected_count = len(collected_coins_total)

# 2. At the apex of his jump, how many coins does Sonic collect on the same
#    vertical level as the highest coin he reaches?
apex_level_count = 0
if coins_collected_at_apex_frame:
    # Find the y-coordinate of the highest coin collected at the apex
    highest_coin_y_at_apex = -1
    for coin in coins_collected_at_apex_frame:
        if coin[1] > highest_coin_y_at_apex:
            highest_coin_y_at_apex = coin[1]
    
    # Count how many coins collected at the apex share this y-level
    for coin in coins_collected_at_apex_frame:
        if coin[1] == highest_coin_y_at_apex:
            apex_level_count += 1

# Print the final answers separated by a space
print(f"{total_collected_count} {apex_level_count}")