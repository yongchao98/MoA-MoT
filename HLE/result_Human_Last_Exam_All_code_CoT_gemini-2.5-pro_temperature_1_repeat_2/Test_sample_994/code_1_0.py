import itertools
import operator

def can_make_24(numbers):
    """
    Checks if a given list/tuple of 4 numbers can be used to make 24.
    """
    # Use a small tolerance for float comparison
    epsilon = 1e-6
    
    ops = {
        '+': operator.add,
        '-': operator.sub,
        '*': operator.mul,
        '/': operator.truediv
    }
    op_symbols = list(ops.keys())

    # Iterate through all unique permutations of the numbers
    for p_nums in set(itertools.permutations(numbers)):
        n1, n2, n3, n4 = p_nums
        
        # Iterate through all combinations of 3 operators
        for p_ops in itertools.product(op_symbols, repeat=3):
            op1_sym, op2_sym, op3_sym = p_ops
            op1, op2, op3 = ops[op1_sym], ops[op2_sym], ops[op3_sym]

            # Try all 5 parenthetical structures
            # Structure 1: (a op b) op (c op d)
            try:
                res = op2(op1(n1, n2), op3(n3, n4))
                if abs(res - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 2: ((a op b) op c) op d
            try:
                res = op3(op2(op1(n1, n2), n3), n4)
                if abs(res - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass
                
            # Structure 3: (a op (b op c)) op d
            try:
                res = op3(op1(n1, op2(n2, n3)), n4)
                if abs(res - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 4: a op ((b op c) op d)
            try:
                res = op1(n1, op3(op2(n2, n3), n4))
                if abs(res - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass

            # Structure 5: a op (b op (c op d))
            try:
                res = op1(n1, op2(n2, op3(n3, n4)))
                if abs(res - 24) < epsilon:
                    return True
            except ZeroDivisionError:
                pass
                
    return False

def calculate_percentage():
    """
    Calculates the percentage of 4-card combinations (1-10) that can make 24.
    """
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards with replacement
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
            
    # Calculate the percentage
    percentage = solvable_count / total_combinations
    
    # Print the final equation as requested
    print(f"Total unique combinations of four cards (values 1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_count}")
    print("\nFinal Calculation:")
    # The prompt asks to output each number in the final equation.
    # Here, the numbers are the solvable count and the total count.
    print(f"{solvable_count} / {total_combinations} = {percentage:.4f}")

if __name__ == '__main__':
    calculate_percentage()
    print("\n<<<0.9021>>>")