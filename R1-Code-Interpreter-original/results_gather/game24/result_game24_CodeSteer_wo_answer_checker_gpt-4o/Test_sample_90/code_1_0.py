from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try all possible ways to insert parentheses
    a, b, c, d = nums
    op1, op2, op3 = ops
    
    # Different ways to group the operations
    expressions = [
        f"(({a}{op1}{b}){op2}{c}){op3}{d}",
        f"({a}{op1}({b}{op2}{c})){op3}{d}",
        f"({a}{op1}{b}){op2}({c}{op3}{d})",
        f"{a}{op1}(({b}{op2}{c}){op3}{d})",
        f"{a}{op1}({b}{op2}({c}{op3}{d}))"
    ]
    
    for expr in expressions:
        try:
            if eval(expr) == 24:
                return expr
        except ZeroDivisionError:
            continue
    return None

def find_solution(numbers):
    operations = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for ops in product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

numbers = [4, 6, 6, 8]
solution = find_solution(numbers)
print(f"<<<{solution}>>>")