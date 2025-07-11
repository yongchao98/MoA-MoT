import itertools
from fractions import Fraction

# Use a global cache (memoization) to store the results of subproblems.
# The key will be a sorted tuple of numbers, and the value will be a set of all possible results.
memo = {}

def get_possible_results(numbers):
    """
    Recursively finds all possible values that can be generated from a list of numbers
    by applying all valid arithmetic operations and groupings.
    """
    # Use a sorted tuple as the key to ensure that (1, 2) and (2, 1) are treated as the same subproblem.
    key = tuple(sorted(numbers))
    if key in memo:
        return memo[key]

    # Base case: If there's only one number, its only possible result is itself.
    if len(numbers) == 1:
        return {numbers[0]}

    results = set()
    
    # Partition the list of numbers into two non-empty sub-lists in every possible way.
    # We iterate up to half the length to avoid duplicate partitions (e.g., {a} & {b,c} is the same as {b,c} & {a}).
    for i in range(1, len(numbers) // 2 + 1):
        for subset1_tuple in itertools.combinations(numbers, i):
            subset1 = list(subset1_tuple)
            
            # Create subset2 from the remaining numbers. This handles duplicates correctly.
            subset2 = list(numbers)
            for num in subset1:
                subset2.remove(num)

            # Recursively find all results for each sub-list.
            res1 = get_possible_results(subset1)
            res2 = get_possible_results(subset2)

            # Combine the results from the two sub-lists using all arithmetic operations.
            for r1 in res1:
                for r2 in res2:
                    results.add(r1 + r2)
                    results.add(r1 * r2)
                    results.add(r1 - r2)
                    results.add(r2 - r1)
                    if r2 != 0:
                        results.add(r1 / r2)
                    if r1 != 0:
                        results.add(r2 / r1)

    # Cache the result for the current set of numbers before returning.
    memo[key] = results
    return results

def calculate_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10)
    that can be used to make 24.
    """
    # Use Fraction for exact arithmetic, with card values from 1 to 10.
    card_values = [Fraction(i) for i in range(1, 11)]

    # Generate all unique sets of 4 cards (combinations with repetition).
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)

    solvable_count = 0
    target = Fraction(24)

    # Check each combination to see if 24 is a possible result.
    for combo in all_combinations:
        possible_values = get_possible_results(list(combo))
        if target in possible_values:
            solvable_count += 1
    
    # Calculate the final percentage.
    percentage = solvable_count / total_combinations

    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"Percentage of solvable combinations: {percentage:.4f}")

# Run the calculation and print the results.
calculate_24_game_percentage()