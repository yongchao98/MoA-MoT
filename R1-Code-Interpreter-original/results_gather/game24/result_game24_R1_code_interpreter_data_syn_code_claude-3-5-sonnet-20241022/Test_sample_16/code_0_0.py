from itertools import permutations
import operator

def get_all_expressions(nums):
    ops = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', operator.mul: '*', operator.truediv: '/'}
    
    def evaluate(nums, ops):
        # Try different parentheses patterns
        try:
            # Pattern 1: ((a op1 b) op2 c) op3 d
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[1](val1, nums[2])
            result1 = ops[2](val2, nums[3])
            if abs(result1 - 24) < 0.0001:
                return f"(({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}{nums[2]}){op_symbols[ops[2]]}{nums[3]}"
            
            # Pattern 2: (a op1 (b op2 c)) op3 d
            val1 = ops[1](nums[1], nums[2])
            val2 = ops[0](nums[0], val1)
            result2 = ops[2](val2, nums[3])
            if abs(result2 - 24) < 0.0001:
                return f"({nums[0]}{op_symbols[ops[0]]}({nums[1]}{op_symbols[ops[1]]}{nums[2]})){op_symbols[ops[2]]}{nums[3]}"
            
            # Pattern 3: a op1 (b op2 (c op3 d))
            val1 = ops[2](nums[2], nums[3])
            val2 = ops[1](nums[1], val1)
            result3 = ops[0](nums[0], val2)
            if abs(result3 - 24) < 0.0001:
                return f"{nums[0]}{op_symbols[ops[0]]}({nums[1]}{op_symbols[ops[1]]}({nums[2]}{op_symbols[ops[2]]}{nums[3]})"
            
        except ZeroDivisionError:
            pass
        return None

    # Try all possible number permutations and operator combinations
    numbers = list(permutations(nums))
    for nums in numbers:
        for op1 in ops:
            for op2 in ops:
                for op3 in ops:
                    result = evaluate(nums, [op1, op2, op3])
                    if result:
                        return result
    return "No solution found"

print(get_all_expressions([6, 7, 12, 12]))