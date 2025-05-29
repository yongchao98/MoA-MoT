from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression with the given numbers and operations
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
                if nums[i + 1] == 0:
                    return None
                result /= nums[i + 1]
        return result
    except ZeroDivisionError:
        return None

def find_expression_to_24(numbers):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            if evaluate_expression(num_perm, ops) == 24:
                # Format the expression as a string
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return f"<<<{expression}>>>"
    return "No solution found"

numbers = [10, 13, 13, 13]
print(find_expression_to_24(numbers))