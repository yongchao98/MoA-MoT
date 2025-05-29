# Let's try different combinations of the numbers to reach 814
from itertools import permutations, product

numbers = [2, 2, 24, 43, 25]
target = 814
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    expression = f"({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]}) {ops[3]} {nums[4]}"
    try:
        if eval(expression) == target:
            return expression
    except ZeroDivisionError:
        pass
    return None

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=4):
        result = evaluate_expression(num_perm, ops)
        if result:
            print(result)
            break