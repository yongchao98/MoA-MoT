from itertools import permutations, product
import operator

def evaluate(nums, ops, paren_pattern):
    if paren_pattern == 0:  # ((a#b)#c)#d
        try:
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[1](val1, nums[2])
            return ops[2](val2, nums[3])
        except:
            return None
    elif paren_pattern == 1:  # (a#b)#(c#d)
        try:
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[2](nums[2], nums[3])
            return ops[1](val1, val2)
        except:
            return None
    elif paren_pattern == 2:  # (a#(b#c))#d
        try:
            val1 = ops[1](nums[1], nums[2])
            val2 = ops[0](nums[0], val1)
            return ops[2](val2, nums[3])
        except:
            return None
    elif paren_pattern == 3:  # a#(b#(c#d))
        try:
            val1 = ops[2](nums[2], nums[3])
            val2 = ops[1](nums[1], val1)
            return ops[0](nums[0], val2)
        except:
            return None
    elif paren_pattern == 4:  # a#((b#c)#d)
        try:
            val1 = ops[1](nums[1], nums[2])
            val2 = ops[2](val1, nums[3])
            return ops[0](nums[0], val2)
        except:
            return None

def find_24(numbers):
    operators = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', 
                 operator.mul: '*', operator.truediv: '/'}
    
    for nums in permutations(numbers):
        for ops in product(operators, repeat=3):
            for pattern in range(5):
                result = evaluate(nums, ops, pattern)
                if result is not None and abs(result - 24) < 1e-10:
                    # Convert to readable format
                    if pattern == 0:
                        return f"(({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}{nums[2]}){op_symbols[ops[2]]}{nums[3]}"
                    elif pattern == 1:
                        return f"({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}({nums[2]}{op_symbols[ops[2]]}{nums[3]})"
                    elif pattern == 2:
                        return f"({nums[0]}{op_symbols[ops[0]]}({nums[1]}{op_symbols[ops[1]]}{nums[2]})){op_symbols[ops[2]]}{nums[3]}"
                    elif pattern == 3:
                        return f"{nums[0]}{op_symbols[ops[0]]}({nums[1]}{op_symbols[ops[1]]}({nums[2]}{op_symbols[ops[2]]}{nums[3]}))"
                    else:
                        return f"{nums[0]}{op_symbols[ops[0]]}(({nums[1]}{op_symbols[ops[1]]}{nums[2]}){op_symbols[ops[2]]}{nums[3]})"
    return "No solution found"

print(find_24([1, 2, 4, 12]))