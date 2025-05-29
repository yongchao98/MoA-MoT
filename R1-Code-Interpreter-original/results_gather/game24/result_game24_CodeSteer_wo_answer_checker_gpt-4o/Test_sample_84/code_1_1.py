from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    try:
        # ((a op1 b) op2 c) op3 d
        if eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}") == 24:
            return f"((({nums[0]} {ops[0]} {nums[1]}) {ops[1]} {nums[2]}) {ops[2]} {nums[3]})"
        # (a op1 (b op2 c)) op3 d
        if eval(f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}") == 24:
            return f"({nums[0]} {ops[0]} ({nums[1]} {ops[1]} {nums[2]})) {ops[2]} {nums[3]}"
        # (a op1 b) op2 (c op3 d)
        if eval(f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})") == 24:
            return f"({nums[0]} {ops[0]} {nums[1]}) {ops[1]} ({nums[2]} {ops[2]} {nums[3]})"
        # a op1 ((b op2 c) op3 d)
        if eval(f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})") == 24:
            return f"{nums[0]} {ops[0]} (({nums[1]} {ops[1]} {nums[2]}) {ops[2]} {nums[3]})"
        # a op1 (b op2 (c op3 d))
        if eval(f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))") == 24:
            return f"{nums[0]} {ops[0]} ({nums[1]} {ops[1]} ({nums[2]} {ops[2]} {nums[3]}))"
    except ZeroDivisionError:
        pass
    return None

def find_expression_for_24(numbers):
    operations = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for ops in product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result} = 24>>>"
    return "No solution found."

# Define the numbers
numbers = [2, 9, 9, 13]

# Find and print the expression
print(find_expression_for_24(numbers))