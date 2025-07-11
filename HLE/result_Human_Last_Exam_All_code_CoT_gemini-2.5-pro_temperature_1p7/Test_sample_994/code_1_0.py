import itertools
from fractions import Fraction

def can_make_24(nums):
    """
    Checks if a list/tuple of 4 numbers can be used to make 24.
    Uses brute force to check all permutations, operator combinations,
    and parenthetical arrangements.
    """
    # Use Python's Fraction for exact arithmetic, avoiding float issues.
    f_nums = [Fraction(n) for n in nums]
    ops = ['+', '-', '*', '/']
    target = 24

    # Using a set gets only unique permutations for lists with repeated numbers
    for p in set(itertools.permutations(f_nums)):
        a, b, c, d = p

        # Iterate through all 4*4*4=64 operator combinations
        for op1, op2, op3 in itertools.product(ops, repeat=3):

            def apply_op(x, y, op):
                """Applies an operator, returning None for division by zero."""
                if op == '+': return x + y
                if op == '-': return x - y
                if op == '*': return x * y
                if op == '/':
                    return x / y if y != 0 else None
                return None

            # Check all 5 parenthetical patterns (Catalan number C_3=5)
            
            # Pattern 1: ((a op b) op c) op d
            res1 = apply_op(a, b, op1)
            if res1 is not None:
                res1 = apply_op(res1, c, op2)
                if res1 is not None:
                    res1 = apply_op(res1, d, op3)
                    if res1 == target: return True

            # Pattern 2: (a op b) op (c op d)
            res2_left = apply_op(a, b, op1)
            res2_right = apply_op(c, d, op3)
            if res2_left is not None and res2_right is not None:
                res2 = apply_op(res2_left, res2_right, op2)
                if res2 == target: return True

            # Pattern 3: (a op (b op c)) op d
            res3 = apply_op(b, c, op2)
            if res3 is not None:
                res3 = apply_op(a, res3, op1)
                if res3 is not None:
                    res3 = apply_op(res3, d, op3)
                    if res3 == target: return True
            
            # Pattern 4: a op ((b op c) op d)
            res4 = apply_op(b, c, op2)
            if res4 is not None:
                res4 = apply_op(res4, d, op3)
                if res4 is not None:
                    res4 = apply_op(a, res4, op1)
                    if res4 == target: return True

            # Pattern 5: a op (b op (c op d))
            res5 = apply_op(c, d, op3)
            if res5 is not None:
                res5 = apply_op(b, res5, op2)
                if res5 is not None:
                    res5 = apply_op(a, res5, op1)
                    if res5 == target: return True

    return False

def main():
    """
    Main function to calculate the percentage of solvable 24-game combinations.
    """
    card_values = range(1, 11)
    
    # Get all unique combinations of 4 cards with values from 1-10
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_combinations = 0
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_combinations += 1
            
    # Calculate and format the percentage
    percentage = solvable_combinations / total_combinations
    
    print(f"Total possible combinations of four card values (1-10): {total_combinations}")
    print(f"Combinations that can be used to make 24: {solvable_combinations}")
    print(f"Percentage of solvable combinations: {percentage:.4f}")

if __name__ == '__main__':
    main()