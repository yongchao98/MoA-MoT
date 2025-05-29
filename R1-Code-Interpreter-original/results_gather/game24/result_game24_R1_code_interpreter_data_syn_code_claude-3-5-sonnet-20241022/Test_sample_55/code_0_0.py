from itertools import permutations
import operator

def evaluate(nums, ops, paren):
    if paren == 0:  # ((a#b)#c)#d
        try:
            r1 = ops[0](nums[0], nums[1])
            r2 = ops[1](r1, nums[2])
            return ops[2](r2, nums[3])
        except:
            return None
    elif paren == 1:  # (a#b)#(c#d)
        try:
            r1 = ops[0](nums[0], nums[1])
            r2 = ops[2](nums[2], nums[3])
            return ops[1](r1, r2)
        except:
            return None
    elif paren == 2:  # (a#(b#c))#d
        try:
            r1 = ops[1](nums[1], nums[2])
            r2 = ops[0](nums[0], r1)
            return ops[2](r2, nums[3])
        except:
            return None

def find_24(numbers):
    operators = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', 
                 operator.mul: '*', operator.truediv: '/'}
    
    for nums in permutations(numbers):
        for op1 in operators:
            for op2 in operators:
                for op3 in operators:
                    ops = [op1, op2, op3]
                    for paren in range(3):
                        result = evaluate(nums, ops, paren)
                        if result is not None and abs(result - 24) < 1e-10:
                            if paren == 0:
                                expr = f"(({nums[0]}{op_symbols[op1]}{nums[1]}){op_symbols[op2]}{nums[2]}){op_symbols[op3]}{nums[3]}"
                            elif paren == 1:
                                expr = f"({nums[0]}{op_symbols[op1]}{nums[1]}){op_symbols[op2]}({nums[2]}{op_symbols[op3]}{nums[3]})"
                            else:
                                expr = f"({nums[0]}{op_symbols[op1]}({nums[1]}{op_symbols[op2]}{nums[2]})){op_symbols[op3]}{nums[3]}"
                            return expr
    return "No solution found"

print(find_24([6, 10, 10, 10]))