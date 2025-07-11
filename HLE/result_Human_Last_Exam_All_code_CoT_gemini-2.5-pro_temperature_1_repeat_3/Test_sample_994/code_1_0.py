import itertools

def can_make_24(nums):
    """
    Recursively checks if a list of numbers can form the number 24 using the
    standard arithmetic operations (+, -, *, /). This function handles all
    permutations of numbers and all groupings of operations (parentheses).
    """
    # Base case: If only one number is left in the list, check if it's 24.
    # We use a small tolerance (epsilon) for floating-point comparisons.
    if len(nums) == 1:
        return abs(nums[0] - 24) < 1e-6

    # Recursive step:
    # Iterate through all unique pairs of numbers in the current list.
    for i in range(len(nums)):
        for j in range(i + 1, len(nums)):
            # Get the pair of numbers.
            a, b = nums[i], nums[j]
            
            # Create a new list containing the remaining numbers.
            remaining_nums = [nums[k] for k in range(len(nums)) if k != i and k != j]
            
            # Apply all possible operations to the pair (a, b) and recurse.
            # The order matters for subtraction and division.
            possible_results = [a + b, a - b, b - a, a * b]
            if b != 0:
                possible_results.append(a / b)
            if a != 0:
                possible_results.append(b / a)
                
            for result in possible_results:
                # If any recursive call finds a solution, propagate True up.
                if can_make_24(remaining_nums + [result]):
                    return True
                    
    # If no combination of pairs and operations yields 24, return False.
    return False

def calculate_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to make 24.
    """
    # The card values are integers from 1 to 10.
    card_values = range(1, 11)

    # Generate all unique combinations of 4 cards, with replacement allowed.
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)

    solvable_count = 0
    # Check each unique combination to see if it can make 24.
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1

    # Calculate the final percentage.
    percentage = solvable_count / total_combinations

    # As requested, output the numbers used in the final percentage calculation.
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can be solved to make 24: {solvable_count}")
    print("The final equation for the percentage is:")
    print(f"{solvable_count} / {total_combinations} = {percentage:.4f}")

# Execute the main function to find and print the result.
calculate_24_game_percentage()