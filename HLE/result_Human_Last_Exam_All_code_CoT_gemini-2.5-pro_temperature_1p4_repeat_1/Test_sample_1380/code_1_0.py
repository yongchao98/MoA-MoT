import random

def simulate_brownian_motion_walk(start_price=100, num_steps=15):
    """
    Simulates and prints a simple 1D random walk, the discrete model
    for Brownian motion, which Bachelier applied to stock prices.

    Args:
        start_price (int): The starting price of the stock.
        num_steps (int): The number of time steps to simulate.
    """
    position = start_price
    print(f"Starting Price: {position}")
    print("-" * 30)

    # The possible steps are -1 (down) or +1 (up)
    possible_steps = [-1, 1]

    for i in range(num_steps):
        # Randomly choose a step
        step = random.choice(possible_steps)
        
        # Store the old position to print the full equation
        old_position = position
        
        # Update the position with the new step
        position += step
        
        # Print each number in the step-by-step calculation
        # The format is: previous_price + (step) = new_price
        print(f"Step {i+1:02d}: {old_position} + ({step}) = {position}")

    print("-" * 30)
    print(f"Final Price after {num_steps} steps: {position}")

# Run the simulation
simulate_brownian_motion_walk()