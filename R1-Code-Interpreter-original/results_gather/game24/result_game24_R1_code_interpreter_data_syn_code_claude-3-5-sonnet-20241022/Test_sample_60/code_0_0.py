from itertools import permutations, product
import math

def evaluate(nums, ops):
    # Try different ways of parentheses placement
    def try_parentheses(numbers, operations):
        try:
            # Case 1: ((a op1 b) op2 c) op3 d
            expr1 = eval(f"({numbers[0]}{operations[0]}{numbers[1]}){operations[1]}{numbers[2]}{operations[2]}{numbers[3]}")
            if math.isclose(expr1, 24, rel_tol=1e-10):
                return f"(({numbers[0]}{operations[0]}{numbers[1]}){operations[1]}{numbers[2]}{operations[2]}{numbers[3]})"
            
            # Case 2: (a op1 (b op2 c)) op3 d
            expr2 = eval(f"({numbers[0]}{operations[0]}({numbers[1]}{operations[1]}{numbers[2]})){operations[2]}{numbers[3]}")
            if math.isclose(expr2, 24, rel_tol=1e-10):
                return f"({numbers[0]}{operations[0]}({numbers[1]}{operations[1]}{numbers[2]})){operations[2]}{numbers[3]}"
            
            # Case 3: a op1 ((b op2 c) op3 d)
            expr3 = eval(f"{numbers[0]}{operations[0]}(({numbers[1]}{operations[1]}{numbers[2]}){operations[2]}{numbers[3]})")
            if math.isclose(expr3, 24, rel_tol=1e-10):
                return f"{numbers[0]}{operations[0]}(({numbers[1]}{operations[1]}{numbers[2]}){operations[2]}{numbers[3]})"
            
            # Case 4: a op1 (b op2 (c op3 d))
            expr4 = eval(f"{numbers[0]}{operations[0]}({numbers[1]}{operations[1]}({numbers[2]}{operations[2]}{numbers[3]})")
            if math.isclose(expr4, 24, rel_tol=1e-10):
                return f"{numbers[0]}{operations[0]}({numbers[1]}{operations[1]}({numbers[2]}{operations[2]}{numbers[3]})"
            
        except:
            return None
        return None

    for ops_perm in product(['+', '-', '*', '/'], repeat=3):
        result = try_parentheses(nums, ops_perm)
        if result:
            return result
    return None

numbers = [2, 8, 8, 13]
for nums in permutations(numbers):
    result = evaluate(nums, None)
    if result:
        print(result)
        break