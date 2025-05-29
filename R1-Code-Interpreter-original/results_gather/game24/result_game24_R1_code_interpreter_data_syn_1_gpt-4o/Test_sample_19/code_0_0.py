from itertools import permutations, product

numbers = [4, 4, 11, 13]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Construct the expression with the given numbers and operations
        expression = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        # Evaluate the expression
        return eval(expression)
    except ZeroDivisionError:
        return None

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        result = evaluate_expression(num_perm, ops)
        if result == target:
            # Construct the expression string
            expression = f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))"
            print(expression)
            break