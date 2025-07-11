# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# List to store the coordinates of collected coins
collected_coins = []

# Main simulation loop
for t in range(time_steps):
    # Update velocity based on acceleration and jump
    vel_x += x_accel * delta_t
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    # Update position based on velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity
    vel_y += g * delta_t
    
    # Ground collision check
    if sonic_pos[1] < 10:
        sonic_pos[1], vel_y = 10, 0

    # A list to hold the coins that were not collected in this frame
    remaining_coins = []
    
    for coin_pos in coins:
        # Define Sonic's bounding box
        sonic_x_min, sonic_x_max = sonic_pos[0], sonic_pos[0] + sonic_size[0]
        sonic_y_min, sonic_y_max = sonic_pos[1], sonic_pos[1] + sonic_size[1]

        # Define the coin's bounding box (size [1, 1])
        coin_x_min, coin_x_max = coin_pos[0], coin_pos[0] + 1
        coin_y_min, coin_y_max = coin_pos[1], coin_pos[1] + 1

        # AABB intersection test
        x_overlap = (sonic_x_min < coin_x_max) and (sonic_x_max > coin_x_min)
        y_overlap = (sonic_y_min < coin_y_max) and (sonic_y_max > coin_y_min)

        if x_overlap and y_overlap:
            # If collision occurs, add to collected list
            collected_coins.append(coin_pos)
        else:
            # If no collision, keep the coin for the next frame
            remaining_coins.append(coin_pos)
            
    # Update the list of available coins for the next frame
    coins = remaining_coins

# --- Part 1: Calculate total coins collected ---
total_collected_count = len(collected_coins)

# --- Part 2: Calculate coins collected at the highest vertical level ---
apex_level_coins = 0
if total_collected_count > 0:
    # Find the maximum y-coordinate among all collected coins
    max_y_collected = max(coin[1] for coin in collected_coins)
    
    # Count how many collected coins are at that y-level
    apex_level_coins = sum(1 for coin in collected_coins if coin[1] == max_y_collected)

# Print the final answers separated by a space
print(f"{total_collected_count} {apex_level_coins}")

<<<21 7>>>