from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    result = nums[0]
    for i in range(3):
        if ops[i] == '+':
            result += nums[i + 1]
        elif ops[i] == '-':
            result -= nums[i + 1]
        elif ops[i] == '*':
            result *= nums[i + 1]
        elif ops[i] == '/':
            if nums[i + 1] == 0 or result % nums[i + 1] != 0:
                return None
            result //= nums[i + 1]
    return result

def find_expression(nums):
    # Try all permutations of numbers and operations
    for num_perm in permutations(nums):
        for ops in product('+-*/', repeat=3):
            if evaluate_expression(num_perm, ops) == 24:
                # Construct the expression string
                expression = f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
                return expression
    return None

numbers = [4, 6, 13, 13]
expression = find_expression(numbers)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")