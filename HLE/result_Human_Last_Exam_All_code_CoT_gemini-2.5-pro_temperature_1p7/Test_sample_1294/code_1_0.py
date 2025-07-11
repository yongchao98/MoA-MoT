import sys

# Here is a minimal Sonic-like game simulation.
# --- Initial Game State and Parameters ---
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
coin_size = [1.0, 1.0]
vel_x, vel_y = 0.0, 0.0
# Generate the grid of coins
coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# --- Simulation & Data Collection ---

# Use a set to store collected coins as tuples (x, y) to ensure uniqueness
collected_coins_total = set()

# Track Sonic's maximum height and position at that moment (apex)
max_sonic_y = -1.0
pos_at_apex = None

# Main simulation loop
for t in range(time_steps):
    # Update velocity and position based on physics
    vel_x += x_accel * delta_t
    # Jump input occurs for a short duration (t=16, 17, 18, 19)
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    vel_y += g * delta_t
    
    # Ground collision check
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0

    # Track the apex of the jump (maximum height)
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]
        # Store a copy of the position at the apex
        pos_at_apex = list(sonic_pos)

    # --- AABB Collision Detection for Total Coins ---
    # Define Sonic's bounding box for the current frame
    sonic_x_min = sonic_pos[0] - sonic_size[0] / 2.0
    sonic_x_max = sonic_pos[0] + sonic_size[0] / 2.0
    sonic_y_min = sonic_pos[1] - sonic_size[1] / 2.0
    sonic_y_max = sonic_pos[1] + sonic_size[1] / 2.0

    # Iterate through all coins to check for new collections
    for coin in coins:
        # Define coin's bounding box
        coin_x_min = coin[0] - coin_size[0] / 2.0
        coin_x_max = coin[0] + coin_size[0] / 2.0
        coin_y_min = coin[1] - coin_size[1] / 2.0
        coin_y_max = coin[1] + coin_size[1] / 2.0

        # AABB intersection test
        if (sonic_x_min < coin_x_max and sonic_x_max > coin_x_min and
            sonic_y_min < coin_y_max and sonic_y_max > coin_y_min):
            collected_coins_total.add(tuple(coin))

# --- Analysis ---

# Question 1: How many coins does Sonic collect in total?
total_collected_count = len(collected_coins_total)

# Question 2: At the apex of his jump, how many coins does Sonic collect 
# on the same vertical level as the highest coin he reaches?
apex_coins_on_highest_level = 0

# Proceed only if any coins were collected and an apex was recorded
if total_collected_count > 0 and pos_at_apex is not None:
    # Find the y-coordinate of the highest coin collected during the entire run
    highest_collected_coin_y = 0.0
    for c in collected_coins_total:
        if c[1] > highest_collected_coin_y:
            highest_collected_coin_y = c[1]

    # Define Sonic's bounding box at the moment of his apex
    apex_sonic_x_min = pos_at_apex[0] - sonic_size[0] / 2.0
    apex_sonic_x_max = pos_at_apex[0] + sonic_size[0] / 2.0
    apex_sonic_y_min = pos_at_apex[1] - sonic_size[1] / 2.0
    apex_sonic_y_max = pos_at_apex[1] + sonic_size[1] / 2.0

    # Count how many coins on that specific y-level he intersects with at the apex
    for coin in coins:
        # We only care about coins on the same level as the highest one collected
        if coin[1] == highest_collected_coin_y:
            coin_x_min = coin[0] - coin_size[0] / 2.0
            coin_x_max = coin[0] + coin_size[0] / 2.0
            coin_y_min = coin[1] - coin_size[1] / 2.0
            coin_y_max = coin[1] + coin_size[1] / 2.0

            # AABB intersection test at the apex
            if (apex_sonic_x_min < coin_x_max and apex_sonic_x_max > coin_x_min and
                apex_sonic_y_min < coin_y_max and apex_sonic_y_max > coin_y_min):
                apex_coins_on_highest_level += 1

# Print the final answers as two numbers separated by a space
print(f"{total_collected_count} {apex_coins_on_highest_level}")
