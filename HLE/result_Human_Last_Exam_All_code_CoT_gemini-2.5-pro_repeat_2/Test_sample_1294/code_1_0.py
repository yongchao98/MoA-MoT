# Here is a minimal Sonic-like game simulation.
# This script is designed to solve a two-part problem based on the simulation:
# 1. Calculate the total number of coins Sonic collects.
# 2. Find the number of coins collected on the same vertical level as the highest coin Sonic reaches.

def solve_sonic_coin_problem():
    """
    Runs the game simulation and calculates the required metrics.
    """
    # --- Simulation Setup ---
    delta_t, time_steps = 0.05, 1000
    g, jump_power, x_accel = -9.8, 75, 5.5

    # Player (Sonic) initial state and properties
    sonic_pos, sonic_size = [5.0, 10.0], [1.0, 2.0]
    vel_x, vel_y = 0.0, 0.0

    # Collectibles (Coins) setup
    # Coins are defined as points, but for collision they have a size of [1,1]
    coins = [[float(x), float(y)] for x in range(15, 46, 2) for y in range(11, 24, 2)]
    coin_size = [1.0, 1.0]

    # --- Tracking and Collision Variables ---
    collected_coin_positions = []

    # Pre-calculate the sum of half-sizes for the AABB collision check.
    # A collision occurs if the distance between centers is less than this sum on both axes.
    x_collision_distance = sonic_size[0] / 2.0 + coin_size[0] / 2.0  # 0.5 + 0.5 = 1.0
    y_collision_distance = sonic_size[1] / 2.0 + coin_size[1] / 2.0  # 1.0 + 0.5 = 1.5

    # --- Main Simulation Loop ---
    for t in range(time_steps):
        # Update velocity based on constant acceleration and jump input
        vel_x += x_accel * delta_t
        # A short burst of upward power simulates a jump
        if 15 < t < 20:
            vel_y += jump_power * delta_t

        # Update position based on current velocity
        sonic_pos[0] += vel_x * delta_t
        sonic_pos[1] += vel_y * delta_t

        # Apply gravity to vertical velocity for the next frame
        vel_y += g * delta_t

        # Check for collision with the ground (at y=10) and reset if necessary
        if sonic_pos[1] < 10.0:
            sonic_pos[1], vel_y = 10.0, 0.0

        # --- AABB Collision Detection with Coins ---
        # Iterate over a copy of the list (coins[:]) because we modify the original list.
        for coin in coins[:]:
            is_colliding_x = abs(sonic_pos[0] - coin[0]) < x_collision_distance
            is_colliding_y = abs(sonic_pos[1] - coin[1]) < y_collision_distance

            if is_colliding_x and is_colliding_y:
                # On collision, store the coin's position and remove it from the simulation.
                collected_coin_positions.append(coin)
                coins.remove(coin)

    # --- Post-Simulation Analysis ---

    # 1. Total number of collected coins.
    total_collected_coins = len(collected_coin_positions)

    # 2. Number of coins on the same vertical level as the highest coin reached.
    coins_on_highest_level = 0
    if total_collected_coins > 0:
        # Find the maximum y-coordinate among all collected coins.
        highest_coin_y = max(c[1] for c in collected_coin_positions)

        # Count how many collected coins share this highest y-coordinate.
        coins_on_highest_level = sum(1 for c in collected_coin_positions if c[1] == highest_coin_y)

    # --- Final Output ---
    # Print the two required numbers separated by a space.
    print(f"{total_collected_coins} {coins_on_highest_level}")

solve_sonic_coin_problem()