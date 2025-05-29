from itertools import permutations
import operator

def get_all_expressions(nums):
    ops = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', operator.mul: '*', operator.truediv: '/'}
    
    def evaluate(nums, ops):
        # Different ways to group operations
        # ((a op1 b) op2 c) op3 d
        try:
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[1](val1, nums[2])
            val3 = ops[2](val2, nums[3])
            if abs(val3 - 24) < 0.0001:
                return f"(({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}{nums[2]}){op_symbols[ops[2]]}{nums[3]}"
        except:
            pass

        # (a op1 b) op2 (c op3 d)
        try:
            val1 = ops[0](nums[0], nums[1])
            val2 = ops[2](nums[2], nums[3])
            val3 = ops[1](val1, val2)
            if abs(val3 - 24) < 0.0001:
                return f"({nums[0]}{op_symbols[ops[0]]}{nums[1]}){op_symbols[ops[1]]}({nums[2]}{op_symbols[ops[2]]}{nums[3]})"
        except:
            pass

        return None

    for num_perm in permutations(nums):
        for op1 in ops:
            for op2 in ops:
                for op3 in ops:
                    result = evaluate(num_perm, [op1, op2, op3])
                    if result:
                        print(result)
                        return

numbers = [10, 13, 13, 13]
get_all_expressions(numbers)