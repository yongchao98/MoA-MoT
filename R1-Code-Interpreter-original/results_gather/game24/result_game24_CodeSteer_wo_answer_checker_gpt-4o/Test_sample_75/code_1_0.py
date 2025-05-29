from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try all possible ways to parenthesize the expression
    a, b, c, d = nums
    op1, op2, op3 = ops
    
    # Different ways to parenthesize
    expressions = [
        f"(({a}{op1}{b}){op2}({c}{op3}{d}))",
        f"(({a}{op1}({b}{op2}{c})){op3}{d})",
        f"({a}{op1}(({b}{op2}{c}){op3}{d}))",
        f"(({a}{op1}{b}){op2}{c}){op3}{d}",
        f"({a}{op1}({b}{op2}({c}{op3}{d})))"
    ]
    
    for expr in expressions:
        try:
            if abs(eval(expr) - 24) < 1e-6:  # Check if the result is 24
                return expr
        except ZeroDivisionError:
            continue
    return None

def find_expression_to_24(numbers):
    for nums in permutations(numbers):
        for ops in product("+-*/", repeat=3):
            expr = evaluate_expression(nums, ops)
            if expr:
                return f"<<<{expr}>>>"
    return "No solution found"

# Example usage
numbers = [2, 8, 12, 12]
print(find_expression_to_24(numbers))