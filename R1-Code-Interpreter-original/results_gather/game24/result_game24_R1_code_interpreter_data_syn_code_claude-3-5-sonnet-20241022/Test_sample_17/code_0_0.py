from itertools import permutations, product
import operator

def evaluate(nums, ops, pattern):
    if pattern == 0:  # ((a#b)#c)#d
        try:
            r1 = ops[0](nums[0], nums[1])
            r2 = ops[1](r1, nums[2])
            return ops[2](r2, nums[3])
        except:
            return None
    elif pattern == 1:  # (a#b)#(c#d)
        try:
            r1 = ops[0](nums[0], nums[1])
            r2 = ops[1](nums[2], nums[3])
            return ops[2](r1, r2)
        except:
            return None
    elif pattern == 2:  # (a#(b#c))#d
        try:
            r1 = ops[1](nums[1], nums[2])
            r2 = ops[0](nums[0], r1)
            return ops[2](r2, nums[3])
        except:
            return None
    return None

def find_24(numbers):
    operators = {
        '+': operator.add,
        '-': operator.sub,
        '*': operator.mul,
        '/': operator.truediv
    }
    
    for nums in permutations(numbers):
        for ops in product([operators[op] for op in '+-*/'], repeat=3):
            for pattern in range(3):
                result = evaluate(nums, ops, pattern)
                if result is not None and abs(result - 24) < 1e-10:
                    op_symbols = {operator.add: '+', operator.sub: '-', 
                                operator.mul: '*', operator.truediv: '/'}
                    if pattern == 0:
                        return f"(({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}{nums[2]}){op_symbols[ops[2]]}{nums[3]}"
                    elif pattern == 1:
                        return f"({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[2]]}({nums[2]}{op_symbols[ops[1]]}{nums[3]})"
                    else:
                        return f"({nums[0]}{op_symbols[ops[0]]}({nums[1]}{op_symbols[ops[1]]}{nums[2]})){op_symbols[ops[2]]}{nums[3]}"
    return "No solution found"

print(find_24([2, 3, 12, 12]))