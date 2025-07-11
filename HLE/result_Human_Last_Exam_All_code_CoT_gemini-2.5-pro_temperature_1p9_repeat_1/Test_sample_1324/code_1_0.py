import math
import numpy as np

def chocolate_game_simulation(initial_chocolates, steps):
    """
    Simulates the chocolate passing game and prints the state at each step.

    Args:
        initial_chocolates (list): A list of even integers for the initial chocolates.
        steps (int): The number of minutes to simulate.
    """
    chocolates = np.array(initial_chocolates, dtype=np.int64)
    n = len(chocolates)

    if n <= 1:
        print("Number of people n must be > 1.")
        return

    # Validate that initial chocolates are all even
    if not all(c % 2 == 0 for c in chocolates):
        print("All initial chocolate counts must be even.")
        return

    print(f"Running simulation for n = {n} people.")
    print(f"Initial state c^0: {chocolates}\n")

    for i in range(steps + 1):
        h = np.max(chocolates)
        l = np.min(chocolates)
        d = h - l

        print(f"--- Minute i = {i} ---")
        print(f"Chocolates c^{i}: {chocolates}")
        print(f"Max h^{i} = {h}")
        print(f"Min l^{i} = {l}")
        print(f"Difference d^{i} = {d}\n")
        
        # Stop if equilibrium is reached
        if d == 0 and i > 0:
            print("Equilibrium reached.")
            break
            
        if i == steps:
            break

        # Store the current state to calculate the next one
        current_chocolates = chocolates.copy()
        for k in range(n):
            # Person k keeps half of their own chocolates and receives half from their left neighbor
            left_neighbor_idx = (k - 1 + n) % n
            
            new_total = (current_chocolates[k] / 2) + (current_chocolates[left_neighbor_idx] / 2)
            
            # If the resulting number is odd, they take an extra chocolate
            if new_total % 2 != 0:
                chocolates[k] = int(new_total + 1)
            else:
                chocolates[k] = int(new_total)

# --- Main execution of the script ---
if __name__ == "__main__":
    # Define an initial state for n=4 people.
    # Choose a state that is not in equilibrium to observe the dynamics.
    initial_setup = [40, 4, 10, 6] 
    
    # Simulate for a number of steps, e.g., 2*n
    simulation_steps = 2 * len(initial_setup)
    
    chocolate_game_simulation(initial_setup, simulation_steps)
