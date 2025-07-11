import itertools
import operator

def can_make_24(numbers):
    """
    Checks if a given tuple of 4 numbers can be used to make 24.
    It exhaustively tries all permutations of numbers, all combinations of
    operators, and all groupings of operations (parentheses).
    """
    ops = {'+': operator.add, '-': operator.sub, '*': operator.mul, '/': operator.truediv}
    epsilon = 1e-6

    # Iterate through all unique permutations of the input numbers
    for p_nums in set(itertools.permutations(numbers)):
        n1, n2, n3, n4 = p_nums
        
        # Iterate through all 4^3 = 64 combinations of 3 operators
        for op_keys in itertools.product(ops.keys(), repeat=3):
            f1, f2, f3 = (ops[k] for k in op_keys)

            # Note: To print the actual equation, you would format a string like:
            # eq_str = f"(({n1} {op_keys[0]} {n2}) {op_keys[1]} {n3}) {op_keys[2]} {n4}"
            # This function returns True on the first solution found for efficiency.

            # Test structure 1: (a op b) op (c op d)
            try:
                if abs(f2(f1(n1, n2), f3(n3, n4)) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Test structure 2: ((a op b) op c) op d
            try:
                if abs(f3(f2(f1(n1, n2), n3), n4) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Test structure 3: (a op (b op c)) op d
            try:
                if abs(f3(f1(n1, f2(n2, n3)), n4) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Test structure 4: a op ((b op c) op d)
            try:
                if abs(f1(n1, f3(f2(n2, n3), n4)) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Test structure 5: a op (b op (c op d))
            try:
                if abs(f1(n1, f2(n2, f3(n3, n4))) - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

    return False

def calculate_solvable_percentage():
    """
    Calculates the percentage of 4-card combinations that can make 24.
    """
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards (values 1-10) with replacement
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
            
    # Calculate and print the result
    percentage = solvable_count / total_combinations
    
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print(f"Percentage of solvable combinations: {percentage:.4f}")

if __name__ == '__main__':
    calculate_solvable_percentage()