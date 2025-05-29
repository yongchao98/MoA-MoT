import itertools

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    try:
        # ((a op1 b) op2 c) op3 d
        if eval(f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}") == 24:
            return f"((({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]})"
        # (a op1 (b op2 c)) op3 d
        if eval(f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}") == 24:
            return f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}"
        # (a op1 b) op2 (c op3 d)
        if eval(f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})") == 24:
            return f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
        # a op1 ((b op2 c) op3 d)
        if eval(f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})") == 24:
            return f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})"
        # a op1 (b op2 (c op3 d))
        if eval(f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))") == 24:
            return f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
    except ZeroDivisionError:
        pass
    return None

def find_expression_to_24(numbers):
    operations = ['+', '-', '*', '/']
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result}>>>"
    return "No solution found."

numbers = [1, 6, 8, 9]
print(find_expression_to_24(numbers))