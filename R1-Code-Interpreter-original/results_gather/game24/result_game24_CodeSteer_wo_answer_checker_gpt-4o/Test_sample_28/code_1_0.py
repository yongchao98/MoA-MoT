from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    try:
        if eval(f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})") == 24:
            return f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        if eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}") == 24:
            return f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
        if eval(f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}") == 24:
            return f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}"
        if eval(f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})") == 24:
            return f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})"
        if eval(f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))") == 24:
            return f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
    except ZeroDivisionError:
        return None
    return None

def find_expression(nums):
    operations = ['+', '-', '*', '/']
    for num_perm in permutations(nums):
        for ops in product(operations, repeat=3):
            expression = evaluate_expression(num_perm, ops)
            if expression:
                return expression
    return "No solution found"

# Input numbers
numbers = [4, 6, 8, 8]

# Find and print the expression
expression = find_expression(numbers)
print(f"<<<{expression}>>>")