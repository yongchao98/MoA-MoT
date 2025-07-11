import math

# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5, 10], [1, 2]
vel_x, vel_y = 0, 0
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]
coin_size = [1, 1]

# Use a set to store collected coins to avoid duplicates.
# Coins are stored as tuples (x, y) to be hashable.
collected_coins = set()

# Main simulation loop
for t in range(time_steps):
    # Update horizontal velocity
    vel_x += x_accel * delta_t
    
    # Apply jump power during a specific time window
    if 15 < t < 20:
        vel_y += jump_power * delta_t
        
    # Update position based on current velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to vertical velocity
    vel_y += g * delta_t
    
    # Check for and resolve collision with the ground
    if sonic_pos[1] < 10:
        sonic_pos[1], vel_y = 10, 0
        
    # Check for coin collection
    sonic_x, sonic_y = sonic_pos
    sonic_w, sonic_h = sonic_size
    
    for coin in coins:
        coin_x, coin_y = coin
        coin_w, coin_h = coin_size
        
        # AABB (Axis-Aligned Bounding Box) intersection test
        is_colliding = (sonic_x < coin_x + coin_w and
                        sonic_x + sonic_w > coin_x and
                        sonic_y < coin_y + coin_h and
                        sonic_y + sonic_h > coin_y)
        
        if is_colliding:
            collected_coins.add(tuple(coin))

# --- Post-simulation analysis ---

# 1. How many coins does Sonic collect in total?
total_collected_count = len(collected_coins)

# 2. How many coins does Sonic collect on the same vertical level as the highest coin he reaches?
coins_at_apex_level_count = 0
if total_collected_count > 0:
    # Find the y-coordinate of the highest collected coin
    highest_coin_y = max(coin[1] for coin in collected_coins)
    
    # Count how many of the collected coins are at that same highest level
    coins_at_apex_level_count = sum(1 for coin in collected_coins if coin[1] == highest_coin_y)

# Print the final answer in the required format
print(f"{total_collected_count} {coins_at_apex_level_count}")