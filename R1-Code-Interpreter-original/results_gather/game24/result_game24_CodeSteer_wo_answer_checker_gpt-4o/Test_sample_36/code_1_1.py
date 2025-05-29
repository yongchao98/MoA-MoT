from itertools import permutations

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    a, b, c, d = nums
    op1, op2, op3 = ops
    
    # Try different parenthesizations
    try:
        if eval(f"({a}{op1}({b}{op2}({c}{op3}{d})))") == 24:
            return f"({a}{op1}({b}{op2}({c}{op3}{d})))"
        if eval(f"(({a}{op1}{b}){op2}({c}{op3}{d}))") == 24:
            return f"(({a}{op1}{b}){op2}({c}{op3}{d}))"
        if eval(f"(({a}{op1}({b}{op2}{c})){op3}{d})") == 24:
            return f"(({a}{op1}({b}{op2}{c})){op3}{d})"
        if eval(f"({a}{op1}(({b}{op2}{c}){op3}{d}))") == 24:
            return f"({a}{op1}(({b}{op2}{c}){op3}{d}))"
        if eval(f"(({a}{op1}{b}){op2}{c}){op3}{d}") == 24:
            return f"(({a}{op1}{b}){op2}{c}){op3}{d}"
    except ZeroDivisionError:
        pass
    return None

def find_expression(nums):
    # Define possible operations
    operations = ['+', '-', '*', '/']
    
    # Generate all permutations of the numbers
    for perm in permutations(nums):
        # Generate all combinations of operations
        for op1 in operations:
            for op2 in operations:
                for op3 in operations:
                    result = evaluate_expression(perm, (op1, op2, op3))
                    if result:
                        return f"<<<{result}>>>"
    return "No solution found"

# Input numbers
numbers = [3, 8, 12, 12]

# Find and print the expression
print(find_expression(numbers))