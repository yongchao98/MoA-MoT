import itertools
import math

def can_make_24(numbers):
    """
    Recursively checks if a list of numbers can be combined to make 24.
    Uses a small tolerance for floating-point comparisons.
    """
    # Base case: if we have one number left, check if it's 24
    if len(numbers) == 1:
        return math.isclose(numbers[0], 24)

    # Recursive step: pick any two distinct numbers, operate on them,
    # and recurse with the new list of numbers.
    # Using permutations is a straightforward way to try both a-b and b-a.
    for a, b in itertools.permutations(numbers, 2):
        remaining_numbers = list(numbers)
        remaining_numbers.remove(a)
        remaining_numbers.remove(b)

        # Try all 4 operations
        # 1. Addition
        if can_make_24([a + b] + remaining_numbers):
            return True
        # 2. Subtraction
        if can_make_24([a - b] + remaining_numbers):
            return True
        # 3. Multiplication
        if can_make_24([a * b] + remaining_numbers):
            return True
        # 4. Division (check for division by zero)
        if b != 0:
            if can_make_24([a / b] + remaining_numbers):
                return True

    return False

def solve_puzzle():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to achieve exactly 24.
    """
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards with replacement
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_combinations = 0
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_combinations += 1
            
    # Calculate the percentage
    percentage = solvable_combinations / total_combinations
    
    # Per the instructions, outputting the numbers in the final equation.
    # The final equation here is the percentage calculation itself.
    print(f"Total unique combinations: {total_combinations}")
    print(f"Combinations that can make 24: {solvable_combinations}")
    print(f"Final Calculation: {solvable_combinations} / {total_combinations} = {percentage:.4f}")

if __name__ == '__main__':
    solve_puzzle()