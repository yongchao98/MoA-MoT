from itertools import permutations, product
import operator

def evaluate(nums, ops, pattern):
    try:
        if pattern == 0:  # ((a#b)#c)#d
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[1](val1, nums[2])
            return ops[2](val2, nums[3])
        elif pattern == 1:  # (a#b)#(c#d)
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[2](nums[2], nums[3])
            return ops[1](val1, val2)
        elif pattern == 2:  # (a#(b#c))#d
            val1 = ops[1](nums[1], nums[2])
            val2 = ops[0](nums[0], val1)
            return ops[2](val2, nums[3])
    except ZeroDivisionError:
        return None

def find_24(numbers):
    operators = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', 
                 operator.mul: '*', operator.truediv: '/'}
    
    for nums in permutations(numbers):
        for ops in product(operators, repeat=3):
            for pattern in range(3):
                result = evaluate(nums, ops, pattern)
                if result is not None and abs(result - 24) < 1e-10:
                    # Convert to readable format
                    if pattern == 0:
                        expr = f"(({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}{nums[2]}){op_symbols[ops[2]]}{nums[3]}"
                    elif pattern == 1:
                        expr = f"({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}({nums[2]}{op_symbols[ops[2]]}{nums[3]})"
                    else:
                        expr = f"({nums[0]}{op_symbols[ops[0]]}({nums[1]}{op_symbols[ops[1]]}{nums[2]})){op_symbols[ops[2]]}{nums[3]}"
                    return expr
    return "No solution found"

print(find_24([2, 2, 2, 11]))