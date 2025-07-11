import random

def find_bachelier_process():
    """
    This function identifies the physical process Louis Bachelier connected
    to financial markets and demonstrates its underlying principle.
    """
    # The name of the physical process Bachelier used as a model.
    process_name = "Brownian motion"

    print(f"The physical process Louis Bachelier discovered was mathematically identical to the spread of stock option probabilities is called: {process_name}.")
    print("\nHe modeled this using a 'random walk', where the next position is the previous position plus a random step.")
    print("Here is a simple simulation of a 1D random walk:\n")

    # Initializing the simulation
    position = 0
    num_steps = 5

    print(f"Initial Position: {position}")

    # Simulating the random walk and printing the "equation" for each step
    for i in range(1, num_steps + 1):
        old_position = position
        # The step can be +1 or -1, chosen randomly
        step = random.choice([-1, 1])
        position = old_position + step
        
        # Printing each number in the equation for the current step
        # Format: New Position = Old Position + (Step)
        print(f"Step {i}: {position} = {old_position} + ({step})")

if __name__ == "__main__":
    find_bachelier_process()