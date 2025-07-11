import math

# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5

# Use floating point numbers for positions to maintain precision during simulation.
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
coin_size = [1.0, 1.0]

# --- Simulation and Analysis ---

# A set to store unique collected coins
collected_coins = set()

# Variables to track the apex of the jump
max_y_pos = -1.0
sonic_pos_at_apex = [0.0, 0.0]

# Pre-calculate collision distances for AABB check (half_width_A + half_width_B)
x_collision_dist = sonic_size[0] / 2 + coin_size[0] / 2
y_collision_dist = sonic_size[1] / 2 + coin_size[1] / 2

# Main simulation loop
for t in range(time_steps):
    # The physics updates are in the same order as the original code snippet.
    vel_x += x_accel * delta_t
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    vel_y += g * delta_t
    
    # Ground check
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0
    
    # After all physics for the frame are resolved, check for state.
    
    # Track the apex by storing the position at the maximum height achieved.
    if sonic_pos[1] > max_y_pos:
        max_y_pos = sonic_pos[1]
        sonic_pos_at_apex = list(sonic_pos)

    # Perform AABB collision check against all coins in every frame.
    for coin in coins:
        if abs(sonic_pos[0] - coin[0]) < x_collision_dist and \
           abs(sonic_pos[1] - coin[1]) < y_collision_dist:
            collected_coins.add(tuple(coin))

# --- Post-simulation Calculations ---

# Part 1: How many coins does Sonic collect in total?
total_collected_count = len(collected_coins)

# Part 2: At the apex, how many coins are collected on the highest collected level?
apex_level_collections = 0
if collected_coins:
    # First, find the y-coordinate of the highest coin Sonic ever collected.
    highest_collected_coin_y = max(c[1] for c in collected_coins)

    # Then, check for collisions AT THE APEX, only with coins on that specific level.
    for coin in coins:
        if coin[1] == highest_collected_coin_y:
            # Perform AABB check using the stored apex position.
            if abs(sonic_pos_at_apex[0] - coin[0]) < x_collision_dist and \
               abs(sonic_pos_at_apex[1] - coin[1]) < y_collision_dist:
                apex_level_collections += 1

# Print the final answers as two numbers separated by a space.
print(f"{total_collected_count} {apex_level_collections}")