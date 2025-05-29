from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    expressions = [
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})",
        f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})",
        f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
    ]
    
    for expr in expressions:
        try:
            if eval(expr) == 24:
                return expr
        except ZeroDivisionError:
            continue
    return None

def find_solution(numbers):
    for nums in permutations(numbers):
        for ops in product("+-*/", repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

numbers = [5, 5, 7, 8]
solution = find_solution(numbers)
print(f"<<<{solution}>>>")