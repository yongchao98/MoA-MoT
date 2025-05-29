from itertools import permutations
import operator

def get_all_expressions(numbers):
    ops = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', operator.mul: '*', operator.truediv: '/'}
    
    def evaluate(nums, ops_list):
        # Try different parentheses arrangements
        try:
            # ((a op1 b) op2 c) op3 d
            val1 = ops_list[0](nums[0], nums[1])
            val2 = ops_list[1](val1, nums[2])
            val3 = ops_list[2](val2, nums[3])
            if abs(val3 - 24) < 0.0001:
                return f"(({nums[0]}{op_symbols[ops_list[0]]}{nums[1]}){op_symbols[ops_list[1]]}{nums[2]}){op_symbols[ops_list[2]]}{nums[3]}"
            
            # (a op1 (b op2 c)) op3 d
            val1 = ops_list[1](nums[1], nums[2])
            val2 = ops_list[0](nums[0], val1)
            val3 = ops_list[2](val2, nums[3])
            if abs(val3 - 24) < 0.0001:
                return f"({nums[0]}{op_symbols[ops_list[0]]}({nums[1]}{op_symbols[ops_list[1]]}{nums[2]})){op_symbols[ops_list[2]]}{nums[3]}"
            
            # a op1 ((b op2 c) op3 d)
            val1 = ops_list[1](nums[1], nums[2])
            val2 = ops_list[2](val1, nums[3])
            val3 = ops_list[0](nums[0], val2)
            if abs(val3 - 24) < 0.0001:
                return f"{nums[0]}{op_symbols[ops_list[0]]}(({nums[1]}{op_symbols[ops_list[1]]}{nums[2]}){op_symbols[ops_list[2]]}{nums[3]})"
            
        except ZeroDivisionError:
            pass
        return None

    for nums in permutations(numbers):
        for op1 in ops:
            for op2 in ops:
                for op3 in ops:
                    result = evaluate(nums, [op1, op2, op3])
                    if result:
                        return result
    return None

numbers = [1, 9, 12, 12]
result = get_all_expressions(numbers)
print(result)