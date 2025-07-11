import itertools
import math

def can_make_24(numbers):
    """
    Recursively checks if a list of numbers can be used to make 24.
    This function tries all combinations of numbers and operations.
    """
    # Base case: if only one number is left, check if it's close to 24.
    # We use math.isclose for safe floating-point comparison.
    if len(numbers) == 1:
        return math.isclose(numbers[0], 24)

    # Recursive step:
    # Iterate through all unique pairs of numbers in the current list.
    for i in range(len(numbers)):
        for j in range(i + 1, len(numbers)):
            a = numbers[i]
            b = numbers[j]
            
            # Create a new list with the remaining numbers.
            remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]

            # Apply all operations to the pair (a, b) and recurse.
            # The result of the operation is added to the list of remaining numbers.
            
            # Case 1: Addition (a + b)
            if can_make_24(remaining + [a + b]):
                return True
            # Case 2: Multiplication (a * b)
            if can_make_24(remaining + [a * b]):
                return True
            # Case 3: Subtraction (a - b and b - a)
            if can_make_24(remaining + [a - b]):
                return True
            if can_make_24(remaining + [b - a]):
                return True
            # Case 4: Division (a / b and b / a), checking for division by zero
            if b != 0 and can_make_24(remaining + [a / b]):
                return True
            if a != 0 and can_make_24(remaining + [b / a]):
                return True

    # If no combination of operations works, return False.
    return False

def calculate_percentage():
    """
    Calculates the percentage of 24-game combinations that are solvable
    for card values from 1 to 10.
    """
    # The possible values on the cards are integers from 1 to 10.
    card_values = range(1, 11)
    
    # Generate all unique combinations of four cards (with replacement allowed).
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    
    total_combinations = len(all_combinations)
    solvable_count = 0
    
    # Loop through each combination and check if it can make 24.
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # Calculate the final percentage.
    percentage = solvable_count / total_combinations
    
    # Print the details of the final calculation as an equation.
    print("Final Calculation:")
    print(f"Percentage of Solvable Combinations = (Number of Solvable Combinations) / (Total Combinations)")
    print(f"Percentage = {solvable_count} / {total_combinations}")
    
    # Print the final result formatted as requested.
    print("\nResult:")
    print(f"The percentage of all possible combinations that can make 24 is: {percentage:.4f}")

# Run the calculation.
if __name__ == "__main__":
    calculate_percentage()