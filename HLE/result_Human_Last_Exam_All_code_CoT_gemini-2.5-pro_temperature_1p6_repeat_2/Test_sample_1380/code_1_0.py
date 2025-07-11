import random

def simulate_bachelier_process(start_price, volatility, num_steps):
    """
    Simulates the random walk process analogous to Brownian motion,
    which Bachelier applied to stock prices.

    Args:
        start_price (float): The initial price of the stock.
        volatility (float): A measure of the size of the random steps (standard deviation).
        num_steps (int): The number of time steps to simulate.
    """
    print(f"Simulating a random walk (Brownian Motion) for a stock price.")
    print(f"Start Price: {start_price:.2f}, Volatility: {volatility}, Steps: {num_steps}\n")

    current_price = start_price
    print(f"Step 0: Start = {current_price:.2f}")

    # The "equation" at each step is: New_Price = Old_Price + Random_Change
    for i in range(1, num_steps + 1):
        # Generate a random change from a normal distribution (mean=0)
        # This is the core of the Wiener process / Brownian motion
        change = random.normalvariate(0, volatility)
        
        previous_price = current_price
        current_price += change
        
        # We must output each number in the step-by-step "equation"
        sign = '+' if change >= 0 else '-'
        print(f"Step {i}: {previous_price:.2f} {sign} {abs(change):.2f} = {current_price:.2f}")
        
    print(f"\nFinal price after {num_steps} steps: {current_price:.2f}")


# Parameters for the simulation
initial_price = 100.0
# Volatility represents the standard deviation of the price changes per step
price_volatility = 2.0
simulation_steps = 10

simulate_bachelier_process(initial_price, price_volatility, simulation_steps)