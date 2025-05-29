import itertools
import operator

def apply_operation(a, b, op):
    """Apply an operation to two numbers."""
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        if b != 0:
            return a / b
    return None

def find_expression(numbers, target=24):
    """Find an expression that evaluates to the target using all numbers."""
    operations = ['+', '-', '*', '/']
    
    # Generate all permutations of the numbers
    for perm in itertools.permutations(numbers):
        # Try all combinations of operations
        for ops in itertools.product(operations, repeat=3):
            # Try different groupings
            # ((a op1 b) op2 c) op3 d
            result1 = apply_operation(perm[0], perm[1], ops[0])
            if result1 is not None:
                result2 = apply_operation(result1, perm[2], ops[1])
                if result2 is not None:
                    result3 = apply_operation(result2, perm[3], ops[2])
                    if result3 == target:
                        return f"<<<(({perm[0]} {ops[0]} {perm[1]}) {ops[1]} {perm[2]}) {ops[2]} {perm[3]}>>>"
            
            # (a op1 (b op2 c)) op3 d
            result1 = apply_operation(perm[1], perm[2], ops[1])
            if result1 is not None:
                result2 = apply_operation(perm[0], result1, ops[0])
                if result2 is not None:
                    result3 = apply_operation(result2, perm[3], ops[2])
                    if result3 == target:
                        return f"<<<({perm[0]} {ops[0]} ({perm[1]} {ops[1]} {perm[2]})) {ops[2]} {perm[3]}>>>"
            
            # a op1 ((b op2 c) op3 d)
            result1 = apply_operation(perm[1], perm[2], ops[1])
            if result1 is not None:
                result2 = apply_operation(result1, perm[3], ops[2])
                if result2 is not None:
                    result3 = apply_operation(perm[0], result2, ops[0])
                    if result3 == target:
                        return f"<<<{perm[0]} {ops[0]} (({perm[1]} {ops[1]} {perm[2]}) {ops[2]} {perm[3]})>>>"
            
            # a op1 (b op2 (c op3 d))
            result1 = apply_operation(perm[2], perm[3], ops[2])
            if result1 is not None:
                result2 = apply_operation(perm[1], result1, ops[1])
                if result2 is not None:
                    result3 = apply_operation(perm[0], result2, ops[0])
                    if result3 == target:
                        return f"<<<{perm[0]} {ops[0]} ({perm[1]} {ops[1]} ({perm[2]} {ops[2]} {perm[3]}))>>>"
            
            # (a op1 b) op2 (c op3 d)
            result1 = apply_operation(perm[0], perm[1], ops[0])
            if result1 is not None:
                result2 = apply_operation(perm[2], perm[3], ops[2])
                if result2 is not None:
                    result3 = apply_operation(result1, result2, ops[1])
                    if result3 == target:
                        return f"<<<({perm[0]} {ops[0]} {perm[1]}) {ops[1]} ({perm[2]} {ops[2]} {perm[3]})>>>"
    
    return "No solution found"

# Test the function with the given input
numbers = [1, 6, 7, 10]
print(find_expression(numbers))