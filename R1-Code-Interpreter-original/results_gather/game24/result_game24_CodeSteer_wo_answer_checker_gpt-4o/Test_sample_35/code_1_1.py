import itertools

def evaluate_expression(nums, ops):
    # Try different ways to group the numbers with parentheses
    expressions = [
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
        f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})",
        f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))",
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
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
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

numbers = [1, 2, 6, 12]
solution = find_solution(numbers)
print(f"<<<{solution}>>>")