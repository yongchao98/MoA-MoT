import itertools

def evaluate_expression(nums, ops):
    # Evaluate the expression with the given numbers and operations
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

def find_expression_for_24(numbers):
    operations = ['+', '-', '*', '/']
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

numbers = [8, 9, 10, 13]
solution = find_expression_for_24(numbers)
print(solution)