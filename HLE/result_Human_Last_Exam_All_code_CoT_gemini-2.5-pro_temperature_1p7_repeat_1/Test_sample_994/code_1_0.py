import itertools
from fractions import Fraction

# A cache (memoization) to store results of subproblems to avoid re-computation.
memo = {}

def can_make_24(numbers):
    """
    Checks if a list of numbers can be used to make 24 using +, -, *, /.
    Uses recursion and memoization for efficiency and Fraction for precision.
    """
    # The key for memoization is the sorted tuple of numbers, ensuring that
    # (1, 2, 3) and (3, 2, 1) are treated as the same subproblem.
    key = tuple(sorted(numbers))
    if key in memo:
        return memo[key]

    # Base case: If only one number is left, a solution is found if it's 24.
    if len(numbers) == 1:
        is_solution = (numbers[0] == 24)
        memo[key] = is_solution
        return is_solution

    # Recursive step: Try all combinations of two numbers from the list.
    for a, b in itertools.combinations(numbers, 2):
        # Create a list of the remaining numbers.
        remaining = list(numbers)
        remaining.remove(a)
        remaining.remove(b)

        # Apply all operations to the pair (a, b) and recurse.
        # Note: a+b and a*b are commutative, so we don't need to check b+a or b*a.
        # Subtraction and division are not, so we check both a-b, b-a and a/b, b/a.
        potential_results = {a + b, a * b, a - b, b - a}
        if b != 0:
            potential_results.add(a / b)
        if a != 0:
            potential_results.add(b / a)

        for result in potential_results:
            # Recursive call with the new set of numbers.
            if can_make_24(remaining + [result]):
                memo[key] = True
                return True

    # If no combination worked after trying all pairs, this subproblem is unsolvable.
    memo[key] = False
    return False

def calculate_solvable_percentage():
    """
    Calculates the percentage of 4-card hands (values 1-10) that can make 24.
    """
    solvable_count = 0
    card_values = range(1, 11)

    # Generate all unique combinations of four cards, allowing repeated values.
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_count = len(all_combinations)

    # Iterate through each unique hand and check if it's solvable.
    for combo in all_combinations:
        # Convert numbers to Fraction for precise arithmetic.
        numbers_as_fractions = [Fraction(n) for n in combo]
        if can_make_24(numbers_as_fractions):
            solvable_count += 1
    
    # Calculate the ratio.
    if total_count > 0:
        ratio = solvable_count / total_count
    else:
        ratio = 0
    
    # Print the result as a final equation, as requested.
    print(f"{solvable_count} / {total_count} = {ratio:.4f}")

if __name__ == '__main__':
    calculate_solvable_percentage()