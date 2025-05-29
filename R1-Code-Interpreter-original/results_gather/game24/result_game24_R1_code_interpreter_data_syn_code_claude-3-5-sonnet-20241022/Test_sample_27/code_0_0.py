from itertools import permutations, product
import math

def evaluate(nums, ops, paren):
    try:
        expr = ''
        if paren == 0:  # ((a#b)#c)#d
            expr = f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
        elif paren == 1:  # (a#(b#c))#d
            expr = f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}"
        elif paren == 2:  # a#((b#c)#d)
            expr = f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})"
        elif paren == 3:  # a#(b#(c#d))
            expr = f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]})"
        
        result = eval(expr)
        if abs(result - 24) < 1e-10:  # Handle floating point comparison
            print(f"{expr} = 24")
            return True
    except:
        return False
    return False

numbers = [11, 13, 13, 13]
operators = ['+', '-', '*', '/']
found = False

for nums in permutations(numbers):
    for ops in product(operators, repeat=3):
        for paren in range(4):
            if evaluate(nums, ops, paren):
                found = True
                break
    if found:
        break