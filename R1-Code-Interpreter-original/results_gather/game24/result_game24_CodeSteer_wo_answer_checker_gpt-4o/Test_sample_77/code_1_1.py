from itertools import permutations

def evaluate_expression(nums, ops):
    # Try different groupings of operations
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

def find_expression_for_24(numbers):
    operations = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for ops in permutations(operations, 3):
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return "No solution found"

numbers = [1, 8, 8, 11]
solution = find_expression_for_24(numbers)
print(f"<<<{solution}>>>")