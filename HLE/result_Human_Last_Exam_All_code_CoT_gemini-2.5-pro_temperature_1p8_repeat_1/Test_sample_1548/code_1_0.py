import math
import random

def h(x, y):
    """
    Approximation of the potential kernel a(x).
    It is positive and harmonic (in the continuous sense) away from the origin
    and grows logarithmically, which are the key properties needed.
    We add a small constant to the argument of the log to avoid issues at norm 1
    and ensure positivity in the simulated neighborhood.
    """
    return math.log(x**2 + y**2 + 1)

def run_simulation(steps, start_pos):
    """
    Simulates the Doob's h-transform of a 2D SRW.

    Args:
        steps (int): The number of steps to simulate.
        start_pos (tuple): The (x, y) starting position.
    """
    x, y = start_pos
    visits_to_positive_x_axis = 0

    print(f"Starting simulation for {steps} steps from position {start_pos}.")
    print("The test transient set is the positive x-axis (y=0, x>0).\n")

    for step in range(steps):
        # Define neighbors
        neighbors = [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]
        weights = []

        # Calculate weights for transition probabilities
        # The probability of moving to a neighbor is proportional to h(neighbor)
        for nx, ny in neighbors:
            if nx == 0 and ny == 0:
                # The h-transformed process never visits the origin
                weights.append(0)
            else:
                weights.append(h(nx, ny))

        if sum(weights) == 0:
            print("Walk is trapped. Cannot move from", (x, y))
            break
        
        # Choose the next position based on the weights
        next_pos = random.choices(neighbors, weights=weights, k=1)[0]
        x, y = next_pos

        # Check if the new position is on the positive x-axis
        if y == 0 and x > 0:
            visits_to_positive_x_axis += 1
            # Uncomment the line below to see when visits happen
            # print(f"Step {step+1}: Visited positive x-axis at {(x, y)}")

    print("Simulation finished.")
    print(f"Final position: {(x, y)}")
    print(f"Final distance from origin: {math.sqrt(x**2 + y**2):.2f}")
    print(f"Total visits to the positive x-axis: {visits_to_positive_x_axis}")


if __name__ == '__main__':
    # Simulation parameters
    N_STEPS = 50000
    # Start away from the origin for the log approximation to be more accurate
    STARTING_POSITION = (10, 10)
    
    run_simulation(N_STEPS, STARTING_POSITION)
