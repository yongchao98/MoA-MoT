from itertools import permutations, product

def apply_operation(a, b, op):
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        if b == 0 or a % b != 0:
            return None
        return a // b

def evaluate_expression(nums, ops):
    # Apply the operations in the order given by ops
    result = apply_operation(nums[0], nums[1], ops[0])
    if result is None:
        return None
    result = apply_operation(result, nums[2], ops[1])
    if result is None:
        return None
    result = apply_operation(result, nums[3], ops[2])
    return result

def find_expression(nums):
    # Try all permutations of numbers and operations
    for num_perm in permutations(nums):
        for ops in product('+-*/', repeat=3):
            if evaluate_expression(num_perm, ops) == 24:
                # Construct the expression string
                expression = f"(({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}{num_perm[2]}){ops[2]}{num_perm[3]}"
                return expression
    return None

numbers = [4, 6, 13, 13]
expression = find_expression(numbers)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")