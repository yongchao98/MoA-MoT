import itertools
from math import isclose

def can_make_24(numbers):
    """
    Checks if a given tuple of 4 numbers can be used to make 24.
    It systematically tries all permutations of the numbers, all combinations
    of operators, and all possible parenthesis arrangements.
    """
    # Use a set to handle duplicate permutations automatically
    for p_nums in set(itertools.permutations(numbers)):
        a, b, c, d = p_nums
        
        # Try all 4^3 = 64 operator combinations
        for ops in itertools.product(['+', '-', '*', '/'], repeat=3):
            op1, op2, op3 = ops

            def do_op(x, op, y):
                """Performs a single operation, returning None for division by zero."""
                if op == '+': return x + y
                if op == '-': return x - y
                if op == '*': return x * y
                if op == '/':
                    return x / y if y != 0 else None
                return None

            # Pattern 1: (a op1 b) op2 (c op3 d)
            res1 = do_op(a, op1, b)
            res2 = do_op(c, op3, d)
            if res1 is not None and res2 is not None:
                final_res = do_op(res1, op2, res2)
                if final_res is not None and isclose(final_res, 24):
                    return True

            # Pattern 2: ((a op1 b) op2 c) op3 d
            res1 = do_op(a, op1, b)
            if res1 is not None:
                res2 = do_op(res1, op2, c)
                if res2 is not None:
                    final_res = do_op(res2, op3, d)
                    if final_res is not None and isclose(final_res, 24):
                        return True
            
            # The remaining 3 patterns are just variations of the order of operations
            # for left-associative operators. Re-shuffling these will be redundant for
            # +, -, *, / but we check them for completeness.

            # Pattern 3: (a op1 (b op2 c)) op3 d
            res1 = do_op(b, op2, c)
            if res1 is not None:
                res2 = do_op(a, op1, res1)
                if res2 is not None:
                    final_res = do_op(res2, op3, d)
                    if final_res is not None and isclose(final_res, 24):
                        return True

            # Pattern 4: a op1 ((b op2 c) op3 d)
            res1 = do_op(b, op2, c)
            if res1 is not None:
                res2 = do_op(res1, op3, d)
                if res2 is not None:
                    final_res = do_op(a, op1, res2)
                    if final_res is not None and isclose(final_res, 24):
                        return True

            # Pattern 5: a op1 (b op2 (c op3 d))
            res1 = do_op(c, op3, d)
            if res1 is not None:
                res2 = do_op(b, op2, res1)
                if res2 is not None:
                    final_res = do_op(a, op1, res2)
                    if final_res is not None and isclose(final_res, 24):
                        return True
                        
    return False

def main():
    """
    Main function to calculate the percentage of solvable 24-game combinations.
    """
    all_card_values = range(1, 11)
    # Generate all unique combinations of 4 cards from values 1-10
    all_combinations = list(itertools.combinations_with_replacement(all_card_values, 4))
    total_combinations_count = len(all_combinations)

    solvable_combinations_count = 0
    # Check each combination
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_combinations_count += 1
            
    # Calculate and print the final result
    percentage_solvable = solvable_combinations_count / total_combinations_count
    
    # Per instructions, output the numbers in the final equation
    print(f"{solvable_combinations_count} / {total_combinations_count} = {percentage_solvable:.4f}")

if __name__ == "__main__":
    main()