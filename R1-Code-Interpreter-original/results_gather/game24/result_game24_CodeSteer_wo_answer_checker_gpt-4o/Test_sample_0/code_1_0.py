from itertools import permutations, product

def evaluate_expression(nums, ops):
    expression = f"({nums[0]}{ops[0]}{nums[1]})"
    result = eval(expression)
    for i in range(2, len(nums)):
        expression = f"({expression}{ops[i-1]}{nums[i]})"
        result = eval(expression)
    return result, expression

def find_expression_to_24(numbers):
    for num_perm in permutations(numbers):
        for ops in product("+-*/", repeat=len(numbers)-1):
            try:
                result, expression = evaluate_expression(num_perm, ops)
                if abs(result - 24) < 1e-6:  # Consider floating-point precision
                    return f"<<<{expression}>>>"
            except ZeroDivisionError:
                continue
    return "No solution found"

numbers = [1, 1, 4, 6]
print(find_expression_to_24(numbers))