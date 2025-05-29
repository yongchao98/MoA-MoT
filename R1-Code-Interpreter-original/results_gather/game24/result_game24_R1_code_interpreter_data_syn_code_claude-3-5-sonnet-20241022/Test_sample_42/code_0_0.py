from itertools import permutations, product
import operator

def evaluate(nums, ops, parens):
    # Convert to list for easier manipulation
    expression = list(nums)
    operations = list(ops)
    
    # Apply operations according to parentheses pattern
    if parens == 0:  # ((a#b)#c)#d
        val = eval(f"{expression[0]}{operations[0]}{expression[1]}")
        val = eval(f"{val}{operations[1]}{expression[2]}")
        val = eval(f"{val}{operations[2]}{expression[3]}")
    elif parens == 1:  # (a#(b#c))#d
        val = eval(f"{expression[1]}{operations[1]}{expression[2]}")
        val = eval(f"{expression[0]}{operations[0]}{val}")
        val = eval(f"{val}{operations[2]}{expression[3]}")
    elif parens == 2:  # a#((b#c)#d)
        val = eval(f"{expression[1]}{operations[1]}{expression[2]}")
        val = eval(f"{val}{operations[2]}{expression[3]}")
        val = eval(f"{expression[0]}{operations[0]}{val}")
    elif parens == 3:  # a#(b#(c#d))
        val = eval(f"{expression[2]}{operations[2]}{expression[3]}")
        val = eval(f"{expression[1]}{operations[1]}{val}")
        val = eval(f"{expression[0]}{operations[0]}{val}")
    return val

def find_24(numbers):
    operators = ['+', '-', '*', '/']
    for nums in permutations(numbers):
        for ops in product(operators, repeat=3):
            for parens in range(4):
                try:
                    result = evaluate(nums, ops, parens)
                    if abs(result - 24) < 0.0001:  # Account for floating point arithmetic
                        if parens == 0:
                            return f"(({nums[0]}{ops[0]}{nums[1]}){ops[1]}{nums[2]}){ops[2]}{nums[3]}"
                        elif parens == 1:
                            return f"({nums[0]}{ops[0]}({nums[1]}{ops[1]}{nums[2]})){ops[2]}{nums[3]}"
                        elif parens == 2:
                            return f"{nums[0]}{ops[0]}(({nums[1]}{ops[1]}{nums[2]}){ops[2]}{nums[3]})"
                        else:
                            return f"{nums[0]}{ops[0]}({nums[1]}{ops[1]}({nums[2]}{ops[2]}{nums[3]}))"
                except:
                    continue
    return "No solution found"

print(find_24([1, 5, 5, 12]))