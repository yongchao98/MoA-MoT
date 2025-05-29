from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    a, b, c, d = nums
    op1, op2, op3 = ops
    
    # Try different parenthesizations
    try:
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
    for nums in permutations(numbers):
        for ops in product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

numbers = [2, 3, 12, 12]
solution = find_expression_for_24(numbers)
print(f"<<<{solution}>>>")