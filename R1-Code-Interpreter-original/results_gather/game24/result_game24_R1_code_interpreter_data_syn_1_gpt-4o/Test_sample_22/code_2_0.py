from itertools import permutations, product

numbers = [1, 6, 6, 11]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Evaluate the expression with different groupings
        expr1 = eval(f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})")
        expr2 = eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}")
        expr3 = eval(f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}")
        expr4 = eval(f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})")
        expr5 = eval(f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))")
        return expr1, expr2, expr3, expr4, expr5
    except ZeroDivisionError:
        return None, None, None, None, None

# Try all permutations of numbers and combinations of operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        results = evaluate_expression(num_perm, ops)
        for i, result in enumerate(results):
            if result == target:
                if i == 0:
                    expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                elif i == 1:
                    expression = f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}"
                elif i == 2:
                    expression = f"({num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}{num_perm[2]})){ops[2]}{num_perm[3]}"
                elif i == 3:
                    expression = f"{num_perm[0]}{ops[0]}(({num_perm[1]}{ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]})"
                elif i == 4:
                    expression = f"{num_perm[0]}{ops[0]}({num_perm[1]}{ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))"
                print(expression)
                break