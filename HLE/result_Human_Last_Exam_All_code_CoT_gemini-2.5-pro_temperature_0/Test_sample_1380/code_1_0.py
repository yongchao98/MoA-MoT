import random

def simulate_bachelier_process(initial_price, num_steps):
    """
    Simulates and prints a simple random walk for a stock price.

    This models the core idea of Bachelier's work: that price changes
    are random and unpredictable, similar to Brownian motion.
    """
    print(f"Simulating a random walk based on Bachelier's model.")
    print(f"Initial Price: {initial_price:.2f}")
    print("-" * 40)

    current_price = float(initial_price)

    for i in range(1, num_steps + 1):
        # Generate a random step, either +1 or -1
        step = random.choice([-1, 1])
        previous_price = current_price
        current_price += step

        # Fulfilling the requirement to show the numbers in the equation
        operator = "+" if step > 0 else "-"
        print(f"Step {i:2}: {previous_price:6.2f} {operator} {abs(step):.2f} = {current_price:6.2f}")

    print("-" * 40)
    print(f"Final Price after {num_steps} steps: {current_price:.2f}")

# --- Parameters for the simulation ---
# You can change these values to see different outcomes
start_price = 100
total_steps = 15

# --- Run the simulation ---
simulate_bachelier_process(start_price, total_steps)