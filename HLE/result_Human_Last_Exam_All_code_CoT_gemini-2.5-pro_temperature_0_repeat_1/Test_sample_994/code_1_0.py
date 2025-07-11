import itertools
import operator

# Use a dictionary for memoization to store results for previously seen number sets
# This significantly speeds up the calculation by avoiding redundant computations.
memo = {}

def get_possible_values(numbers):
    """
    Recursively calculates all possible values that can be formed from a list of numbers
    using the four basic arithmetic operations.
    """
    # Use a sorted tuple as a key for the memoization cache because the order of numbers
    # in a sub-problem does not matter for the set of possible results.
    numbers_tuple = tuple(sorted(numbers))
    if numbers_tuple in memo:
        return memo[numbers_tuple]

    # Base case: If there's only one number, that's the only possible value.
    if len(numbers) == 1:
        return {numbers[0]}

    results = set()
    # Partition the set of numbers into two non-empty subsets.
    # We iterate from 1 to half the length of the list to generate unique partitions.
    for i in range(1, len(numbers) // 2 + 1):
        for part1_indices in itertools.combinations(range(len(numbers)), i):
            part1_nums = [numbers[j] for j in part1_indices]
            part2_nums = [numbers[j] for j in range(len(numbers)) if j not in part1_indices]

            # Recursively find all possible values for each partition.
            values1 = get_possible_values(part1_nums)
            values2 = get_possible_values(part2_nums)

            # Combine the results from the two partitions with all operators.
            for v1 in values1:
                for v2 in values2:
                    results.add(v1 + v2)
                    results.add(v1 - v2)
                    results.add(v2 - v1)
                    results.add(v1 * v2)
                    if v2 != 0:
                        results.add(v1 / v2)
                    if v1 != 0:
                        results.add(v2 / v1)
    
    # Cache the result before returning.
    memo[numbers_tuple] = results
    return results

def can_make_24(numbers):
    """
    Checks if the number 24 can be made from a given list of four numbers.
    """
    # Clear the memoization cache for each new top-level combination.
    memo.clear()
    all_values = get_possible_values(numbers)
    # Check if 24 is in the results, using a small tolerance for floating-point comparisons.
    for val in all_values:
        if abs(val - 24) < 1e-6:
            return True
    return False

def calculate_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """
    card_values = range(1, 11)
    # Generate all unique combinations of 4 cards with replacement.
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    # Calculate the percentage and format the output.
    percentage = solvable_count / total_combinations
    
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"The final calculation is: {solvable_count} / {total_combinations} = {percentage:.4f}")

if __name__ == '__main__':
    calculate_24_game_percentage()