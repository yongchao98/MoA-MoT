from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different ways to parenthesize the expression
    a, b, c, d = nums
    op1, op2, op3 = ops
    
    # Evaluate all possible parenthesizations
    expressions = [
        f"(({a}{op1}{b}){op2}({c}{op3}{d}))",
        f"(({a}{op1}({b}{op2}{c})){op3}{d})",
        f"({a}{op1}(({b}{op2}{c}){op3}{d}))",
        f"(({a}{op1}{b}){op2}{c}){op3}{d}",
        f"({a}{op1}({b}{op2}({c}{op3}{d})))"
    ]
    
    for expr in expressions:
        try:
            if abs(eval(expr) - 24) < 1e-9:  # Check if the result is 24
                return expr
        except ZeroDivisionError:
            continue
    return None

def find_equation_to_24(numbers):
    # Generate all permutations of numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product("+-*/", repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result}>>>"
    return "No solution found"

# Test the function with the given input
numbers = [4, 5, 7, 7]
print(find_equation_to_24(numbers))