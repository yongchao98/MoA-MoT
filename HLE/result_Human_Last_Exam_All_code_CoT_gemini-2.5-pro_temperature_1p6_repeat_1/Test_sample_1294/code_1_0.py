import sys
# Here is a minimal Sonic-like game simulation.
# The following code implements the plan to complete the simulation and find the answers.

# Initialize all simulation parameters as specified in the problem
delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
sonic_pos, sonic_size = [5, 10], [1, 2]
vel_x, vel_y = 0, 0

# Generate the initial list of all coins. Coin size is assumed to be [1,1].
coins = [[x, y] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# We use a set to store the coordinates of collected coins.
# A set is efficient and automatically prevents duplicate entries.
collected_coins = set()

# Main simulation loop that runs for each time step.
for t in range(time_steps):
    # Update Sonic's horizontal velocity due to constant acceleration.
    vel_x += x_accel * delta_t
    
    # Apply a vertical velocity boost if Sonic is jumping (within the specified time window).
    if 15 < t < 20:
        vel_y += jump_power * delta_t
    
    # Update Sonic's position using the current velocities.
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t
    
    # Apply gravity to the vertical velocity.
    vel_y += g * delta_t
    
    # Check for collision with the ground (at y=10) and reset position/velocity if so.
    if sonic_pos[1] < 10:
        sonic_pos[1], vel_y = 10, 0
        
    # --- Collision Detection Stage ---
    
    # Define Sonic's bounding box for the current frame.
    sonic_left = sonic_pos[0]
    sonic_right = sonic_pos[0] + sonic_size[0]
    sonic_bottom = sonic_pos[1]
    sonic_top = sonic_pos[1] + sonic_size[1]
    
    # A temporary list to store coins that were NOT collected in this frame.
    remaining_coins = []
    
    # Iterate through all currently available coins to check for collisions.
    for coin_pos in coins:
        coin_x, coin_y = coin_pos[0], coin_pos[1]
        
        # Define the coin's bounding box (size is [1, 1]).
        coin_left = coin_x
        coin_right = coin_x + 1
        coin_bottom = coin_y
        coin_top = coin_y + 1
        
        # Perform the AABB intersection test. A collision occurs if the boxes overlap on both axes.
        is_colliding = (sonic_left < coin_right and sonic_right > coin_left and
                        sonic_bottom < coin_top and sonic_top > coin_bottom)
        
        if is_colliding:
            # If a collision occurs, add the coin's position (as a tuple) to our set of collected coins.
            collected_coins.add(tuple(coin_pos))
        else:
            # If there's no collision, this coin remains in play for the next frame.
            remaining_coins.append(coin_pos)
            
    # Update the main list of coins, effectively removing the ones collected in this frame.
    coins = remaining_coins

# --- Analysis Stage ---

# Part 1: How many coins does Sonic collect in total?
# The answer is the total number of unique coins in our 'collected_coins' set.
total_collected_coins = len(collected_coins)

# Part 2: At the apex of his jump, how many coins does Sonic collect on the same
# vertical level as the highest coin he reaches?
# This is interpreted as: "Of all the coins collected, what is the highest y-level,
# and how many collected coins are at that level?"
coins_on_highest_level = 0
if total_collected_coins > 0:
    # First, find the maximum y-coordinate among all the collected coins.
    highest_coin_y = max(coin[1] for coin in collected_coins)
    
    # Next, count how many of the collected coins share this highest y-coordinate.
    coins_on_highest_level = sum(1 for coin in collected_coins if coin[1] == highest_coin_y)

# Final Output: Print the two calculated numbers, separated by a space.
print(total_collected_coins, coins_on_highest_level)
<<<51 6>>>