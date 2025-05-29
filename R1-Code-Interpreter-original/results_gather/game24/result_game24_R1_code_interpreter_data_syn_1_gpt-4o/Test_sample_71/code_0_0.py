import itertools

def evaluate_expression(nums, ops):
    # Evaluate the expression with the given numbers and operations
    try:
        result = eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))")
        return result
    except ZeroDivisionError:
        return None

def find_expression_to_24(numbers):
    # All possible operations
    operations = ['+', '-', '*', '/']
    
    # Generate all permutations of numbers and operations
    for num_perm in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            if evaluate_expression(num_perm, ops) == 24:
                return f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]}))"
    return None

numbers = [8, 9, 11, 12]
expression = find_expression_to_24(numbers)
print(expression)