import itertools
from math import isclose

def solve_24_game():
    """
    Calculates the percentage of 4-card combinations (values 1-10) that can make 24.
    """

    # Helper function to apply an operator safely, checking for division by zero.
    def apply_op(a, b, op):
        if op == '+': return a + b
        if op == '-': return a - b
        if op == '*': return a * b
        if op == '/':
            if b == 0:
                # Return a value that will never equal 24 to signify failure.
                return float('inf')
            return a / b
        return float('inf')

    # Main solver for a single combination of 4 numbers.
    def can_make_24(numbers):
        """
        Checks if a given list/tuple of 4 numbers can make 24 using +, -, *, /.
        """
        # Iterate through all unique permutations of the numbers.
        # Using set() handles combinations with duplicate numbers correctly and efficiently.
        for p in set(itertools.permutations(numbers)):
            a, b, c, d = p
            
            # Iterate through all 4^3 = 64 operator combinations.
            for op1, op2, op3 in itertools.product(['+', '-', '*', '/'], repeat=3):
                
                # Check all 5 parenthesis patterns for 4 numbers.
                
                # Pattern 1: (a op b) op (c op d)
                res1 = apply_op(a, b, op1)
                res2 = apply_op(c, d, op2)
                if isclose(apply_op(res1, res2, op3), 24):
                    return True
                
                # Pattern 2: ((a op b) op c) op d
                res1 = apply_op(a, b, op1)
                res2 = apply_op(res1, c, op2)
                if isclose(apply_op(res2, d, op3), 24):
                    return True
                
                # Pattern 3: a op (b op (c op d))
                res1 = apply_op(c, d, op3)
                res2 = apply_op(b, res1, op2)
                if isclose(apply_op(a, res2, op1), 24):
                    return True
                    
                # Pattern 4: a op ((b op c) op d)
                res1 = apply_op(b, c, op2)
                res2 = apply_op(res1, d, op3)
                if isclose(apply_op(a, res2, op1), 24):
                    return True

                # Pattern 5: (a op (b op c)) op d
                res1 = apply_op(b, c, op2)
                res2 = apply_op(a, res1, op1)
                if isclose(apply_op(res2, d, op3), 24):
                    return True

        return False

    # Main execution logic
    card_values = range(1, 11)
    
    # Generate all unique combinations of 4 cards with values from 1-10.
    all_combinations = list(itertools.combinations_with_replacement(card_values, 4))
    total_combinations = len(all_combinations)
    
    solvable_count = 0
    for combo in all_combinations:
        if can_make_24(combo):
            solvable_count += 1
            
    # Calculate the final proportion.
    result = solvable_count / total_combinations
    
    # Print the numbers used in the final calculation.
    print(f"Calculation: {solvable_count} / {total_combinations}")
    print(f"Result: {result:.4f}")

if __name__ == '__main__':
    solve_24_game()
