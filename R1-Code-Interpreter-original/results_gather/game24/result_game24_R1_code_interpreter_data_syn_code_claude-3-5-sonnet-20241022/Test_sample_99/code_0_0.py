from itertools import permutations
import operator

def evaluate(nums, ops):
    # Helper function to evaluate expression with different parentheses positions
    def eval_expr(expr):
        try:
            return eval(expr)
        except:
            return None

    n = nums
    o = ops
    
    # Try different parentheses positions
    expressions = [
        f"({n[0]}{o[0]}{n[1]}){o[1]}({n[2]}{o[2]}{n[3]})",
        f"(({n[0]}{o[0]}{n[1]}){o[1]}{n[2]}){o[2]}{n[3]}",
        f"{n[0]}{o[0]}({n[1]}{o[1]}({n[2]}{o[2]}{n[3]}))",
        f"{n[0]}{o[0]}(({n[1]}{o[1]}{n[2]}){o[2]}{n[3]})"
    ]
    
    for expr in expressions:
        result = eval_expr(expr)
        if result is not None and abs(result - 24) < 1e-10:
            return expr
    return None

numbers = [3, 3, 9, 9]
operators = ['+', '-', '*', '/']
found = False

for nums in permutations(numbers):
    for op1 in operators:
        for op2 in operators:
            for op3 in operators:
                expr = evaluate(nums, [op1, op2, op3])
                if expr is not None:
                    print(f"{expr} = 24")
                    found = True
                    break
            if found:
                break
        if found:
            break
    if found:
        break