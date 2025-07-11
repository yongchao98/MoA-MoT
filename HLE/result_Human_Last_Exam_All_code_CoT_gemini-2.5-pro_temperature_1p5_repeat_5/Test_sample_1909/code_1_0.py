import random

def run_simulation(n):
    """
    Runs a single simulation of the crossing-out process for a row of n integers.

    Args:
        n (int): The initial number of integers in the row.

    Returns:
        int: The number of integers remaining at the end of the process.
    """
    if n == 0:
        return 0
    
    # Use a list of booleans to represent the presence of numbers 1 to n.
    # We use size n+1 to have a 1-based indexing for convenience.
    is_present = [True] * (n + 1)

    while True:
        # Find all available consecutive pairs
        consecutive_pairs = []
        for i in range(1, n):
            if is_present[i] and is_present[i+1]:
                consecutive_pairs.append(i)
        
        # If no pairs are left, the process stops
        if not consecutive_pairs:
            break
            
        # Randomly select a pair to remove
        # We stored the first number of the pair, e.g., i for the pair (i, i+1)
        chosen_i = random.choice(consecutive_pairs)
        
        # Remove the chosen pair
        is_present[chosen_i] = False
        is_present[chosen_i+1] = False
        
    # Count the number of remaining integers
    remaining_count = sum(is_present[1:])
    return remaining_count

def main():
    """
    Main function to run multiple simulations and compute the final ratio.
    """
    try:
        n = 1000  # Number of integers in the row (for approximating infinity)
        num_trials = 2000 # Number of simulations to run for a stable average
    except NameError:
        n = 1000
        num_trials = 2000

    total_remaining_count = 0
    for _ in range(num_trials):
        total_remaining_count += run_simulation(n)
        
    avg_remaining_count = total_remaining_count / num_trials
    ratio = avg_remaining_count / n
    
    # Final output as an equation
    print(f"Simulation parameters: n = {n}, trials = {num_trials}")
    print("Expected value of (Number of remaining numbers / n)")
    print(f"= E[R_n] / n")
    print(f"~= ({avg_remaining_count:.4f}) / {n}")
    print(f"~= {ratio:.4f}")

if __name__ == "__main__":
    main()