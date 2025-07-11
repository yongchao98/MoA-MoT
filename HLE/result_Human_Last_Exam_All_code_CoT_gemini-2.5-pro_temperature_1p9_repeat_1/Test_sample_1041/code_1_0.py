import numpy as np
import numba

# --- Configuration ---
START_POS = (0, 300)
TARGET_SET = {(0, 0), (2, 0)}
R_MAX = 1000
R_MAX_SQ = R_MAX**2
N_WALKS = 1_000_000 # Number of simulations

@numba.jit(nopython=True)
def run_simulation(n_walks, start_pos_x, start_pos_y, r_max_sq, target_1_x, target_1_y, target_2_x, target_2_y):
    """
    Runs the Monte Carlo simulation for the random walk.
    This function is accelerated with Numba.
    """
    hits = 0
    # Define the four possible steps in a 2D simple random walk
    steps = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])

    for i in range(n_walks):
        pos_x, pos_y = start_pos_x, start_pos_y

        # We can set a maximum number of steps to prevent any theoretical infinite loops,
        # although escape or hit is guaranteed for a 2D random walk.
        # A large enough number, e.g., 10,000,000, is safe.
        for _ in range(10_000_000):
            # Check for a hit
            if (pos_x == target_1_x and pos_y == target_1_y) or \
               (pos_x == target_2_x and pos_y == target_2_y):
                hits += 1
                break

            # Check for leaving the disk
            if pos_x**2 + pos_y**2 > r_max_sq:
                break
            
            # Take a random step
            step_idx = np.random.randint(0, 4)
            pos_x += steps[step_idx, 0]
            pos_y += steps[step_idx, 1]

    return hits

def main():
    """
    Main function to run the simulation and print the results.
    """
    print("Running Monte Carlo simulation...")
    print("This may take a moment depending on your computer's speed.")
    print("Please ensure you have 'numba' installed for optimal performance (`pip install numba`).")

    # Extract target coordinates for numba compatibility
    target_list = list(TARGET_SET)
    target_1_x, target_1_y = target_list[0]
    target_2_x, target_2_y = target_list[1]
    
    hits = run_simulation(N_WALKS, START_POS[0], START_POS[1], R_MAX_SQ, target_1_x, target_1_y, target_2_x, target_2_y)
    
    probability = hits / N_WALKS
    
    # --- Final Equation Output ---
    # The probability is calculated as the ratio of successful walks to the total number of walks.
    print("\n--- Calculation ---")
    print(f"Number of hits = {hits}")
    print(f"Total number of walks = {N_WALKS}")
    print(f"Probability = Number of hits / Total number of walks")
    print(f"Probability = {hits} / {N_WALKS}")
    
    # Output the final probability with three significant digits
    print("\n--- Result ---")
    print(f"The estimated probability is: {probability:.3g}")


if __name__ == "__main__":
    main()
