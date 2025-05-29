import itertools

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the given numbers and operations
    try:
        result = nums[0]
        for i in range(3):
            if ops[i] == '+':
                result += nums[i + 1]
            elif ops[i] == '-':
                result -= nums[i + 1]
            elif ops[i] == '*':
                result *= nums[i + 1]
            elif ops[i] == '/':
                result /= nums[i + 1]
        return result
    except ZeroDivisionError:
        return None

def find_expression(nums):
    operations = ['+', '-', '*', '/']
    for num_perm in itertools.permutations(nums):
        for ops in itertools.product(operations, repeat=3):
            if evaluate_expression(num_perm, ops) == 24:
                # Format the expression
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return f"<<<{expression}>>>"
    return "No solution found"

nums = [2, 2, 2, 11]
print(find_expression(nums))