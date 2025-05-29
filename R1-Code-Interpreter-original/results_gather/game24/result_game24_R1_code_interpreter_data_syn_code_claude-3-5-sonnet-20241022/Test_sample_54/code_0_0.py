from itertools import permutations
import operator

def evaluate(nums, ops):
    # Helper function to evaluate expression with given numbers and operators
    def calc(a, b, op):
        if op == '+': return a + b
        if op == '-': return a - b
        if op == '*': return a * b
        if op == '/' and b != 0: return a / b
        return None

    # Try different parentheses positions
    # ((a op b) op c) op d
    try:
        res1 = calc(nums[0], nums[1], ops[0])
        if res1 is not None:
            res2 = calc(res1, nums[2], ops[1])
            if res2 is not None:
                res3 = calc(res2, nums[3], ops[2])
                if res3 is not None and abs(res3 - 24) < 0.0001:
                    return f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
    except:
        pass

    # (a op b) op (c op d)
    try:
        res1 = calc(nums[0], nums[1], ops[0])
        res2 = calc(nums[2], nums[3], ops[1])
        if res1 is not None and res2 is not None:
            res3 = calc(res1, res2, ops[2])
            if res3 is not None and abs(res3 - 24) < 0.0001:
                return f"({nums[0]}{ops[0]}{nums[1]}){ops[2]}({nums[2]}{ops[1]}{nums[3]})"
    except:
        pass

    return None

numbers = [3, 7, 9, 9]
operators = ['+', '-', '*', '/']
found = False

for nums in permutations(numbers):
    for op1 in operators:
        for op2 in operators:
            for op3 in operators:
                result = evaluate(nums, [op1, op2, op3])
                if result:
                    print(result + " = 24")
                    found = True
                    break
            if found: break
        if found: break
    if found: break