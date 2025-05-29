import itertools

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
                result /= nums[i + 1]
        return result
    except ZeroDivisionError:
        return None

def find_expression_to_24(numbers):
    # All possible operations
    operations = ['+', '-', '*', '/']
    
    # Generate all permutations of numbers and operations
    for num_perm in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            # Evaluate the expression
            result = evaluate_expression(num_perm, ops)
            if result == 24:
                # Format the expression as a string
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return expression
    return None

numbers = [3, 8, 11, 12]
expression = find_expression_to_24(numbers)
print(expression)