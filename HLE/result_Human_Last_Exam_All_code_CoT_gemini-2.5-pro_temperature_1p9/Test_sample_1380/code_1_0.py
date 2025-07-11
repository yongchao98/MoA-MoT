import random

def simulate_bachelier_walk(start_price, num_days):
    """
    Simulates a simple random walk similar to the process Bachelier described.
    This process is a discrete analogue of Brownian Motion.

    Args:
      start_price (float): The initial price of the stock.
      num_days (int): The number of days to simulate.
    """
    price = start_price
    print(f"Starting Price: {price}")
    print("-" * 30)

    # In this simplified model, we'll represent the process as:
    # Price(t+1) = Price(t) + Change
    # Where 'Change' is a random variable. Here we use +1 or -1.

    for day in range(1, num_days + 1):
        # Generate a random change, either +1 or -1
        change = random.choice([-1, 1])

        # The 'equation' for this step is: New Price = Old Price + Change
        # We print each number in this step's equation.
        old_price = price
        price += change
        
        # Format the change as a signed number for the equation printout
        change_str = f"+ {abs(change)}" if change > 0 else f"- {abs(change)}"
        
        print(f"Day {day:2}: {price:4.0f} = {old_price:4.0f} {change_str}")

# --- Simulation Parameters ---
initial_price = 100
simulation_days = 20

simulate_bachelier_walk(initial_price, simulation_days)