import random

def simulate_random_walk(start_price, num_steps, volatility):
    """
    Simulates a simple random walk similar to the process Bachelier modeled.

    Args:
        start_price (float): The initial price of the asset.
        num_steps (int): The number of time steps to simulate.
        volatility (float): The size of the random up/down move at each step.
    """
    price_path = [start_price]
    current_price = start_price

    print(f"Simulating a random walk for {num_steps} steps, starting at ${start_price:.2f}.")
    print("This process is mathematically identical to Brownian motion.")
    print("-" * 50)
    print("The simulated price path is:")

    # We start with the first price already in the list
    # so we generate the next num_steps prices
    for i in range(num_steps):
        # Generate a random move, either up or down
        move = random.choice([-volatility, volatility])
        current_price += move
        # Ensure price doesn't go below zero
        current_price = max(0, current_price)
        price_path.append(current_price)

    # The prompt asks to "output each number in the final equation".
    # We interpret this as printing the final sequence of prices.
    # We will format them to look like an equation showing the path.
    path_string = " -> ".join([f"{price:.2f}" for price in price_path])
    print(path_string)

# --- Parameters ---
initial_price = 100
simulation_steps = 15
step_size = 2 # Represents the volatility

simulate_random_walk(initial_price, simulation_steps, step_size)