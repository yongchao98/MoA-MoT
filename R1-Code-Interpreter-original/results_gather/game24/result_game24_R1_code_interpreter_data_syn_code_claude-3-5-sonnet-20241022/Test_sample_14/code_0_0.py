from itertools import permutations
import operator

def evaluate(nums, ops):
    # Try different parentheses combinations
    try:
        # ((a op1 b) op2 c) op3 d
        val1 = ops[0](nums[0], nums[1])
        val2 = ops[1](val1, nums[2])
        val3 = ops[2](val2, nums[3])
        if abs(val3 - 24) < 0.0001:
            return True, f"(({nums[0]}{ops[0].__name__}{nums[1]}){ops[1].__name__}{nums[2]}){ops[2].__name__}{nums[3]}"
        
        # (a op1 (b op2 c)) op3 d
        val1 = ops[1](nums[1], nums[2])
        val2 = ops[0](nums[0], val1)
        val3 = ops[2](val2, nums[3])
        if abs(val3 - 24) < 0.0001:
            return True, f"({nums[0]}{ops[0].__name__}({nums[1]}{ops[1].__name__}{nums[2]})){ops[2].__name__}{nums[3]}"
        
        # a op1 ((b op2 c) op3 d)
        val1 = ops[1](nums[1], nums[2])
        val2 = ops[2](val1, nums[3])
        val3 = ops[0](nums[0], val2)
        if abs(val3 - 24) < 0.0001:
            return True, f"{nums[0]}{ops[0].__name__}(({nums[1]}{ops[1].__name__}{nums[2]}){ops[2].__name__}{nums[3]})"
    except:
        pass
    return False, ""

numbers = [7, 7, 11, 12]
operations = [operator.add, operator.sub, operator.mul, operator.truediv]
ops_symbols = {operator.add: '+', operator.sub: '-', operator.mul: '*', operator.truediv: '/'}

for nums in permutations(numbers):
    for op1 in operations:
        for op2 in operations:
            for op3 in operations:
                found, expr = evaluate(nums, [op1, op2, op3])
                if found:
                    # Convert operator names to symbols
                    for op, symbol in ops_symbols.items():
                        expr = expr.replace(op.__name__, symbol)
                    print(expr)
                    exit()