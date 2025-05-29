from itertools import permutations, product

def evaluate_expression(nums, ops):
    """Evaluate the expression given numbers and operations."""
    a, b, c, d = nums
    op1, op2, op3 = ops
    try:
        # Try different parenthesizations
        if eval(f"({a}{op1}{b}){op2}({c}{op3}{d})") == 24:
            return f"(({a}{op1}{b}){op2}({c}{op3}{d}))"
        if eval(f"(({a}{op1}{b}){op2}{c}){op3}{d}") == 24:
            return f"(({a}{op1}{b}){op2}{c}){op3}{d}"
        if eval(f"({a}{op1}({b}{op2}{c})){op3}{d}") == 24:
            return f"({a}{op1}({b}{op2}{c})){op3}{d}"
        if eval(f"{a}{op1}(({b}{op2}{c}){op3}{d})") == 24:
            return f"{a}{op1}(({b}{op2}{c}){op3}{d})"
        if eval(f"{a}{op1}({b}{op2}({c}{op3}{d}))") == 24:
            return f"{a}{op1}({b}{op2}({c}{op3}{d}))"
    except ZeroDivisionError:
        pass
    return None

def find_expression_to_24(numbers):
    """Find an expression that evaluates to 24 using all numbers."""
    operations = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for ops in product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result}>>>"
    return "No solution found."

# Input numbers
numbers = [3, 4, 11, 12]

# Find and print the expression
print(find_expression_to_24(numbers))