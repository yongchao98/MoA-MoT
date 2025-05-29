from itertools import permutations
import operator

def get_operations():
    return [
        (operator.add, '+'),
        (operator.sub, '-'),
        (operator.mul, '*'),
        (operator.truediv, '/')
    ]

def evaluate(nums, ops, pattern):
    try:
        if pattern == 0:  # ((a op1 b) op2 c) op3 d
            temp = ops[0][0](nums[0], nums[1])
            temp = ops[1][0](temp, nums[2])
            return ops[2][0](temp, nums[3])
        elif pattern == 1:  # (a op1 (b op2 c)) op3 d
            temp = ops[1][0](nums[1], nums[2])
            temp = ops[0][0](nums[0], temp)
            return ops[2][0](temp, nums[3])
        elif pattern == 2:  # a op1 (b op2 (c op3 d))
            temp = ops[2][0](nums[2], nums[3])
            temp = ops[1][0](nums[1], temp)
            return ops[0][0](nums[0], temp)
        elif pattern == 3:  # a op1 ((b op2 c) op3 d)
            temp = ops[1][0](nums[1], nums[2])
            temp = ops[2][0](temp, nums[3])
            return ops[0][0](nums[0], temp)
    except ZeroDivisionError:
        return None
    except:
        return None

numbers = [2, 8, 12, 12]
operations = get_operations()
target = 24
found = False

for nums in permutations(numbers):
    for op1 in operations:
        for op2 in operations:
            for op3 in operations:
                ops = [op1, op2, op3]
                for pattern in range(4):
                    result = evaluate(nums, ops, pattern)
                    if result is not None and abs(result - target) < 0.0001:
                        expression = ""
                        if pattern == 0:
                            expression = f"(({nums[0]} {op1[1]} {nums[1]}) {op2[1]} {nums[2]}) {op3[1]} {nums[3]}"
                        elif pattern == 1:
                            expression = f"({nums[0]} {op1[1]} ({nums[1]} {op2[1]} {nums[2]})) {op3[1]} {nums[3]}"
                        elif pattern == 2:
                            expression = f"{nums[0]} {op1[1]} ({nums[1]} {op2[1]} ({nums[2]} {op3[1]} {nums[3]}))"
                        elif pattern == 3:
                            expression = f"{nums[0]} {op1[1]} (({nums[1]} {op2[1]} {nums[2]}) {op3[1]} {nums[3]})"
                        print(f"{expression} = 24")
                        found = True
                        break
                if found:
                    break
            if found:
                break
        if found:
            break
    if found:
        break