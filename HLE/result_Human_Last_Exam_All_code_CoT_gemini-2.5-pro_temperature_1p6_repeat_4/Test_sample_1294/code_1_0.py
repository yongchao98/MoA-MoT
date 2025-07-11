import sys

# Here is a minimal Sonic-like game simulation.
# The script will simulate the game and then perform analysis to answer the questions.

# --- Step 1: Initialization ---
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5

# Use floating-point numbers for position and velocity for accuracy
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0

# Generate coins as a set of tuples for efficient lookup
coins = {(x, y) for x in range(15, 46, 2) for y in range(11, 24, 2)}

# --- Step 2: Prepare for Data Collection ---
# Use a set to store unique collected coins for the total count
collected_coins_total = set()

# Variables to track the jump apex
max_sonic_y = -1.0
history = []

# --- Step 3: Main Simulation Loop ---
for t in range(time_steps):
    # Update velocity based on acceleration and jump inputs
    vel_x += x_accel * delta_t
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    # Update position based on velocity at the start of the frame
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to velocity for the next frame
    vel_y += g * delta_t
    
    # Ground check: if Sonic is below the floor, reset position and vertical velocity
    if sonic_pos[1] < 10.0:
        sonic_pos[1], vel_y = 10.0, 0.0

    # Track the highest point of the jump (the apex)
    if sonic_pos[1] > max_sonic_y:
        max_sonic_y = sonic_pos[1]

    # --- Step 4: Collision Detection and Data Logging ---
    s_x, s_y = sonic_pos[0], sonic_pos[1]
    s_w, s_h = sonic_size[0], sonic_size[1]
    
    coins_collected_this_frame = set()
    
    # Check for collisions with all coins
    for coin_x, coin_y in coins:
        c_x, c_y = float(coin_x), float(coin_y)
        c_w, c_h = 1.0, 1.0  # Coin size is [1, 1]

        # AABB intersection test
        if (s_x < c_x + c_w and
            s_x + s_w > c_x and
            s_y < c_y + c_h and
            s_y + s_h > c_y):
            
            collected_coin = (coin_x, coin_y)
            coins_collected_this_frame.add(collected_coin)

    # Add the coins collected in this frame to the total set
    collected_coins_total.update(coins_collected_this_frame)
    
    # Store frame data for later analysis of the apex collection
    history.append({
        'pos': list(sonic_pos), # Store a copy of the position
        'collected': coins_collected_this_frame
    })

# --- Step 5: Final Analysis ---

# Question 1: How many coins does Sonic collect in total?
total_collected_count = len(collected_coins_total)

# Question 2: At the apex of his jump, how many coins does Sonic collect 
# on the same vertical level as the highest coin he reaches?

apex_level_count = 0
if collected_coins_total:
    # Find the y-coordinate of the highest coin Sonic ever collected
    highest_coin_y = max(c[1] for c in collected_coins_total)
    
    # Identify all unique coins collected during the apex frame(s)
    coins_collected_at_apex = set()
    for frame_data in history:
        # Check if this frame was at the apex height (use a small tolerance for float comparison)
        if abs(frame_data['pos'][1] - max_sonic_y) < 1e-9:
            coins_collected_at_apex.update(frame_data['collected'])

    # Count how many of those collected at the apex have the required y-coordinate
    for coin in coins_collected_at_apex:
        if coin[1] == highest_coin_y:
            apex_level_count += 1

# --- Step 6: Output ---
# Print the two final numbers separated by a space.
print(f"{total_collected_count} {apex_level_count}")