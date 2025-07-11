# Here is a minimal Sonic-like game simulation.

def aabb_intersect(pos1, size1, pos2, size2):
    """
    Checks if two Axis-Aligned Bounding Boxes intersect.
    pos is [x, y], size is [width, height].
    Returns True if they intersect, False otherwise.
    """
    x1, y1 = pos1
    w1, h1 = size1
    x2, y2 = pos2
    w2, h2 = size2

    # Check for no overlap on X axis
    if x1 + w1 < x2 or x2 + w2 < x1:
        return False
    # Check for no overlap on Y axis
    if y1 + h1 < y2 or y2 + h2 < y1:
        return False
    # If there's no gap on either axis, they must overlap.
    return True

# --- Simulation Parameters and Initial State ---
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5, 10], [1, 2]
vel_x, vel_y = 0, 0
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]
coin_size = [1, 1]  # As per the problem description

# --- Analysis Variables ---
collected_coins = set()  # Using a set to store (x, y) tuples of unique collected coins
max_sonic_y = -float('inf')
apex_sonic_pos = None
highest_collected_coin_y = -1

# --- Main Simulation Loop ---
for t in range(time_steps):
    # Update horizontal velocity and apply jump impulse
    vel_x += x_accel * delta_t
    if 15 < t < 20:  # Jump on time steps 16, 17, 18, 19
        vel_y += jump_power * delta_t
    
    # Update position based on velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to vertical velocity
    vel_y += g * delta_t
    
    # Ground collision check
    if sonic_pos[1] < 10:
        sonic_pos[1], vel_y = 10, 0

    # --- In-loop analysis ---
    # Check for the apex of the jump
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]
        apex_sonic_pos = list(sonic_pos) # Store a copy of Sonic's position at his apex

    # Check for coin collection
    for coin_pos in coins:
        if aabb_intersect(sonic_pos, sonic_size, coin_pos, coin_size):
            coin_tuple = tuple(coin_pos)
            if coin_tuple not in collected_coins:
                collected_coins.add(coin_tuple)
                # Update the highest Y coordinate of any collected coin
                if coin_pos[1] > highest_collected_coin_y:
                    highest_collected_coin_y = coin_pos[1]

# --- Post-simulation Calculation ---

# Task 1: Total number of unique coins collected
total_coins_collected = len(collected_coins)

# Task 2: At the apex of the jump, count coins collected on the highest collected coin's level
coins_at_apex_on_level = 0
if apex_sonic_pos is not None and highest_collected_coin_y != -1:
    for coin_pos in coins:
        # Check if the coin is on the correct vertical level
        if coin_pos[1] == highest_collected_coin_y:
            # Check if Sonic collides with this coin at his apex position
            if aabb_intersect(apex_sonic_pos, sonic_size, coin_pos, coin_size):
                coins_at_apex_on_level += 1

# Final output
print(f"{total_coins_collected} {coins_at_apex_on_level}")
