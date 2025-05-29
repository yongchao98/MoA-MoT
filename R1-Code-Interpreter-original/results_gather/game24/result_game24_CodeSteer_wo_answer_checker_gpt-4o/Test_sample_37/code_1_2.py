import itertools

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    try:
        # ((a op1 b) op2 c) op3 d
        expr1 = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
        if eval(expr1) == 24:
            return expr1
        # (a op1 (b op2 c)) op3 d
        expr2 = f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}"
        if eval(expr2) == 24:
            return expr2
        # (a op1 b) op2 (c op3 d)
        expr3 = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
        if eval(expr3) == 24:
            return expr3
        # a op1 ((b op2 c) op3 d)
        expr4 = f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})"
        if eval(expr4) == 24:
            return expr4
        # a op1 (b op2 (c op3 d))
        expr5 = f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
        if eval(expr5) == 24:
            return expr5
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