import itertools

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    expressions = [
        f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})",
        f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}",
        f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}",
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

def find_expression_to_24(numbers):
    operations = ['+', '-', '*', '/']
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            expr = evaluate_expression(nums, ops)
            if expr:
                return expr
    return "No solution found"

numbers = [1, 4, 6, 8]
solution = find_expression_to_24(numbers)
print(f"<<<{solution}>>>")