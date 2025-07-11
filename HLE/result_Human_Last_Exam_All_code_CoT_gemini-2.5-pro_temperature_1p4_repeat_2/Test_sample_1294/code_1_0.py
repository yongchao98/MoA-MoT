# Here is a minimal Sonic-like game simulation.
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5, 10], [1, 2]
vel_x, vel_y = 0, 0
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# Use a set to store the coordinates of collected coins to avoid duplicates.
collected_coins = set()
coin_size = [1, 1]

# The main simulation loop runs for the given number of time steps.
for t in range(time_steps):
    # Update Sonic's horizontal velocity due to acceleration.
    vel_x += x_accel * delta_t
    
    # Apply a vertical impulse for the jump during a specific time window.
    if 15 < t < 20: # This corresponds to t = 16, 17, 18, 19
        vel_y += jump_power * delta_t
    
    # Update Sonic's position using the current velocity.
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to the vertical velocity for the next frame.
    vel_y += g * delta_t
    
    # Handle collision with the ground. If Sonic is below y=10, reset his position and vertical velocity.
    if sonic_pos[1] < 10:
        sonic_pos[1] = 10
        vel_y = 0

    # Define Sonic's bounding box for the current frame.
    sonic_min_x = sonic_pos[0]
    sonic_max_x = sonic_pos[0] + sonic_size[0]
    sonic_min_y = sonic_pos[1]
    sonic_max_y = sonic_pos[1] + sonic_size[1]
    
    # Iterate through all coins to check for collisions.
    for coin_pos in coins:
        # Define the coin's bounding box.
        coin_min_x = coin_pos[0]
        coin_max_x = coin_pos[0] + coin_size[0]
        coin_min_y = coin_pos[1]
        coin_max_y = coin_pos[1] + coin_size[1]
        
        # AABB (Axis-Aligned Bounding Box) intersection test.
        is_colliding = (sonic_min_x < coin_max_x and sonic_max_x > coin_min_x and
                        sonic_min_y < coin_max_y and sonic_max_y > coin_min_y)
        
        if is_colliding:
            # If a collision occurs, add the coin's position to our set of collected coins.
            collected_coins.add(tuple(coin_pos))

# --- Analysis ---

# 1. How many coins does Sonic collect in total?
# This is simply the number of unique coins in our collected_coins set.
total_collected = len(collected_coins)

# 2. At the apex of his jump, how many coins does Sonic collect 
#    on the same vertical level as the highest coin he reaches?
# This is interpreted as: among all the coins collected during the simulation,
# find the maximum y-coordinate, and then count how many collected coins share that y-coordinate.
count_at_highest_level = 0
if total_collected > 0:
    # Find the y-coordinate of the highest coin collected.
    max_collected_y = max(coin[1] for coin in collected_coins)
    
    # Count how many of the collected coins are at that same highest level.
    count_at_highest_level = sum(1 for coin in collected_coins if coin[1] == max_collected_y)

# Print the two answers separated by a space.
print(f"{total_collected} {count_at_highest_level}")