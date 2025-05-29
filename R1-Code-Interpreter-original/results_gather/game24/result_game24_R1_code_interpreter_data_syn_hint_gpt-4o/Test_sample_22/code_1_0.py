from itertools import permutations, product

numbers = [1, 6, 6, 11]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        # Evaluate the expression in the order of operations
        result = nums[0]
        for i in range(3):
            if ops[i] == '+':
                result += nums[i + 1]
            elif ops[i] == '-':
                result -= nums[i + 1]
            elif ops[i] == '*':
                result *= nums[i + 1]
            elif ops[i] == '/':
                if nums[i + 1] == 0:
                    return None
                result /= nums[i + 1]
        return result
    except ZeroDivisionError:
        return None

for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        if evaluate_expression(num_perm, ops) == target:
            expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
            print(expression)
            break