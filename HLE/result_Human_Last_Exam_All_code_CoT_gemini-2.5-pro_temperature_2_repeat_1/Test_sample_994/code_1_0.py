import itertools

# A memoization cache to store results for sets of numbers that have already been computed.
memo = {}

def is_solvable(numbers):
    """
    Recursively checks if a list of numbers can be used to form the number 24
    using addition, subtraction, multiplication, and division.
    """
    # Use a sorted tuple as the key for memoization, as the order of numbers
    # in the list doesn't affect its solvability.
    key = tuple(sorted(numbers))
    if key in memo:
        return memo[key]

    # Base case: If only one number is left, check if it's 24.
    # We use a small tolerance (epsilon) for floating-point comparisons.
    if len(numbers) == 1:
        is_24 = abs(numbers[0] - 24) < 1e-9
        memo[key] = is_24
        return is_24

    # Recursive step: Pick any two numbers from the list and apply an operation.
    # Using set(itertools.combinations) gets all unique pairs.
    for a, b in set(itertools.combinations(numbers, 2)):
        
        # Create a new list containing the numbers that were not picked.
        remaining = list(numbers)
        remaining.remove(a)
        remaining.remove(b)

        # Try all operations on the pair (a, b) and recurse.
        # Note: Subtraction and division are not commutative, so we must
        # test both 'a op b' and 'b op a'.
        
        # Test: a + b
        if is_solvable(remaining + [a + b]):
            memo[key] = True
            return True
            
        # Test: a * b
        if is_solvable(remaining + [a * b]):
            memo[key] = True
            return True
            
        # Test: a - b and b - a
        if is_solvable(remaining + [a - b]) or is_solvable(remaining + [b - a]):
            memo[key] = True
            return True

        # Test: a / b and b / a (checking for division by zero)
        if b != 0 and is_solvable(remaining + [a / b]):
            memo[key] = True
            return True
        if a != 0 and is_solvable(remaining + [b / a]):
            memo[key] = True
            return True

    # If no operations lead to a solution, this set is unsolvable.
    memo[key] = False
    return False

def calculate_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to make 24.
    """
    # 1. Generate all unique combinations of four card values (1-10).
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_count = len(all_combinations)

    # 2. Count how many of these combinations are solvable.
    solvable_count = 0
    for combo in all_combinations:
        # Pass the combination as a list to the solver function.
        if is_solvable(list(combo)):
            solvable_count += 1
            
    # 3. Calculate the final percentage.
    percentage = solvable_count / total_count

    print(f"Total unique combinations of four cards (values 1-10): {total_count}")
    print(f"Combinations that can make 24: {solvable_count}")
    print("Final equation for the percentage:")
    # We output each number in the final equation as requested.
    print(f"{solvable_count} / {total_count} = {percentage:.4f}")

if __name__ == '__main__':
    calculate_24_game_percentage()