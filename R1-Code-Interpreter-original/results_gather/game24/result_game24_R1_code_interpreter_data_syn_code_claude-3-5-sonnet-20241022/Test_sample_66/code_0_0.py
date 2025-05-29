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
            r2 = ops[1](nums[2], nums[3])
            return ops[2](r1, r2)
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
    operators = {
        '+': operator.add,
        '-': operator.sub,
        '*': operator.mul,
        '/': operator.truediv
    }
    
    for nums in permutations(numbers):
        for ops in product(operators.values(), repeat=3):
            for paren in range(3):
                result = evaluate(nums, ops, paren)
                if result is not None and abs(result - 24) < 1e-10:
                    op_symbols = [op.__name__.replace('truediv', '/').replace('sub', '-').replace('add', '+').replace('mul', '*') for op in ops]
                    if paren == 0:
                        return f"(({nums[0]}{op_symbols[0]}{nums[1]}){op_symbols[1]}{nums[2]}){op_symbols[2]}{nums[3]}"
                    elif paren == 1:
                        return f"({nums[0]}{op_symbols[0]}{nums[1]}){op_symbols[2]}({nums[2]}{op_symbols[1]}{nums[3]})"
                    else:
                        return f"({nums[0]}{op_symbols[0]}({nums[1]}{op_symbols[1]}{nums[2]})){op_symbols[2]}{nums[3]}"
    return "No solution found"

print(find_24([5, 7, 10, 12]))