import itertools

numbers = [1, 8, 10, 13]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        result = eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))")
        return result
    except ZeroDivisionError:
        return None

for num_perm in itertools.permutations(numbers):
    for ops in itertools.product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))"
            print(expression)
            break