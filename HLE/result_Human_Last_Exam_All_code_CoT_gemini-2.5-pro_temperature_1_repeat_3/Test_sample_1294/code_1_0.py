import math

# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
# Use floating point numbers for position for better precision
sonic_pos, sonic_size = [5.0, 10.0], [1, 2]
vel_x, vel_y = 0.0, 0.0
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# Variables to track collected coins
collected_coin_indices = set()
collected_coins_positions = []
sonic_w, sonic_h = sonic_size[0], sonic_size[1]
coin_w, coin_h = 1, 1

for t in range(time_steps):
    # Update horizontal velocity and apply jump force
    vel_x += x_accel * delta_t
    if 15 < t < 20:
        vel_y += jump_power * delta_t

    # Update position based on velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to vertical velocity
    vel_y += g * delta_t
    
    # Check for ground collision
    if sonic_pos[1] < 10:
        sonic_pos[1], vel_y = 10, 0
        
    # Check for coin collection
    sonic_x, sonic_y = sonic_pos
    for i, coin in enumerate(coins):
        # Skip coins that have already been collected
        if i in collected_coin_indices:
            continue
        
        coin_x, coin_y = coin
        
        # AABB (Axis-Aligned Bounding Box) intersection test
        is_colliding = (sonic_x < coin_x + coin_w and
                        sonic_x + sonic_w > coin_x and
                        sonic_y < coin_y + coin_h and
                        sonic_y + sonic_h > coin_y)
        
        if is_colliding:
            collected_coin_indices.add(i)
            collected_coins_positions.append(coin)

# --- Part 1: Calculate the total number of collected coins ---
total_collected = len(collected_coin_indices)

# --- Part 2: Calculate coins collected at the highest vertical level ---
coins_at_highest_level = 0
if collected_coins_positions:
    # Find the maximum y-coordinate among all collected coins
    max_y = -1
    for coin in collected_coins_positions:
        if coin[1] > max_y:
            max_y = coin[1]
            
    # Count how many collected coins are at that maximum y-level
    count = 0
    for coin in collected_coins_positions:
        if coin[1] == max_y:
            count += 1
    coins_at_highest_level = count

# Print the final answers
print(f"{total_collected} {coins_at_highest_level}")