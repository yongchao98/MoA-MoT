from itertools import permutations, product
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
    elif paren == 3:  # a#(b#(c#d))
        try:
            r1 = ops[2](nums[2], nums[3])
            r2 = ops[1](nums[1], r1)
            return ops[0](nums[0], r2)
        except:
            return None
    elif paren == 4:  # a#((b#c)#d)
        try:
            r1 = ops[1](nums[1], nums[2])
            r2 = ops[2](r1, nums[3])
            return ops[0](nums[0], r2)
        except:
            return None

def find_24(numbers):
    ops = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', 
                 operator.mul: '*', operator.truediv: '/'}
    
    for nums in permutations(numbers):
        for operators in product(ops, repeat=3):
            for paren in range(5):
                result = evaluate(nums, operators, paren)
                if result is not None and abs(result - 24) < 1e-10:
                    # Convert to readable format
                    if paren == 0:
                        return f"(({nums[0]}{op_symbols[operators[0]]}{nums[1]}){op_symbols[operators[1]]}{nums[2]}){op_symbols[operators[2]]}{nums[3]}"
                    elif paren == 1:
                        return f"({nums[0]}{op_symbols[operators[0]]}{nums[1]}){op_symbols[operators[1]]}({nums[2]}{op_symbols[operators[2]]}{nums[3]})"
                    elif paren == 2:
                        return f"({nums[0]}{op_symbols[operators[0]]}({nums[1]}{op_symbols[operators[1]]}{nums[2]})){op_symbols[operators[2]]}{nums[3]}"
                    elif paren == 3:
                        return f"{nums[0]}{op_symbols[operators[0]]}({nums[1]}{op_symbols[operators[1]]}({nums[2]}{op_symbols[operators[2]]}{nums[3]}))"
                    else:
                        return f"{nums[0]}{op_symbols[operators[0]]}(({nums[1]}{op_symbols[operators[1]]}{nums[2]}){op_symbols[operators[2]]}{nums[3]})"
    return "No solution found"

print(find_24([3, 5, 5, 6]))