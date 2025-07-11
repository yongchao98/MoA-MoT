# Here is a minimal Sonic-like game simulation.

delta_t, time_steps = 0.05, 1000
g, jump_power, x_accel = -9.8, 75, 5.5
# Use floating point numbers for physics calculations to maintain precision
sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
vel_x, vel_y = 0.0, 0.0
coin_size = [1.0, 1.0] # Coins are size [1,1]
# Initialize coins with float coordinates
coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]

# This list will store the coordinates of coins as they are collected
collected_coins = []

# Main simulation loop
for t in range(time_steps):
    # Update horizontal velocity with acceleration
    vel_x += x_accel * delta_t

    # Apply jump power for t = 16, 17, 18, 19, as per the condition '15 < t < 20'
    if 15 < t < 20:
        vel_y += jump_power * delta_t

    # Update Sonic's position based on his velocity
    sonic_pos[0] += vel_x * delta_t
    sonic_pos[1] += vel_y * delta_t

    # Apply gravity to vertical velocity
    vel_y += g * delta_t

    # Check for collision with the ground (at y=10)
    if sonic_pos[1] < 10.0:
        sonic_pos[1] = 10.0
        vel_y = 0.0

    # This list will hold the coins that are NOT collected in the current frame
    remaining_coins = []
    # Check for collisions with each coin
    for coin in coins:
        # Define Sonic's and the coin's bounding boxes for the AABB test
        sx, sy = sonic_pos
        sw, sh = sonic_size
        cx, cy = coin
        cw, ch = coin_size

        # AABB intersection test:
        # True if the bounding boxes overlap on both the X and Y axes
        is_colliding = (sx < cx + cw and
                        sx + sw > cx and
                        sy < cy + ch and
                        sy + sh > cy)

        if is_colliding:
            # If a collision occurs, add the coin to our collected list
            collected_coins.append(coin)
        else:
            # If no collision, the coin remains for the next frame
            remaining_coins.append(coin)

    # Update the list of available coins, effectively removing the ones just collected
    coins = remaining_coins

# --- Analysis after the simulation ---

# 1. How many coins does Sonic collect in total?
total_collected = len(collected_coins)

# 2. How many coins does Sonic collect on the same vertical level as the highest coin he reaches?
apex_level_count = 0
if total_collected > 0:
    # Find the maximum y-coordinate among all the collected coins
    highest_coin_y = max(c[1] for c in collected_coins)

    # Count how many of the collected coins are on that same highest y-level
    apex_level_count = sum(1 for c in collected_coins if c[1] == highest_coin_y)

# Print the two answers separated by a space
print(f"{total_collected} {apex_level_count}")