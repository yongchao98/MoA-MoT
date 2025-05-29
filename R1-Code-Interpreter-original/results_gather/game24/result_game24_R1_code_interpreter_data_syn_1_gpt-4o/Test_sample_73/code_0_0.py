from itertools import permutations, product

numbers = [2, 6, 6, 8]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Construct the expression with the given numbers and operations
        expression = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        # Evaluate the expression
        result = eval(expression)
        # Check if the result is equal to the target
        if result == target:
            return expression
    except ZeroDivisionError:
        return None
    return None

# Try all permutations of numbers and combinations of operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        expr = evaluate_expression(num_perm, ops)
        if expr:
            print(expr)
            break