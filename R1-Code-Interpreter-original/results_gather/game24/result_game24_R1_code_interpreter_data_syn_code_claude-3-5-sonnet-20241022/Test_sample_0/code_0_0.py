from itertools import permutations
import operator

def generate_expressions(numbers):
    ops = [operator.add, operator.sub, operator.mul, operator.truediv]
    op_symbols = {operator.add: '+', operator.sub: '-', 
                  operator.mul: '*', operator.truediv: '/'}
    
    def evaluate(nums, ops_list):
        # Different parentheses patterns
        try:
            # ((a#b)#c)#d
            val1 = ops_list[2](ops_list[1](ops_list[0](nums[0], nums[1]), nums[2]), nums[3])
            if abs(val1 - 24) < 0.0001:
                return f"(({nums[0]}{op_symbols[ops_list[0]]}{nums[1]}){op_symbols[ops_list[1]]}{nums[2]}){op_symbols[ops_list[2]]}{nums[3]}"
            
            # (a#(b#c))#d
            val2 = ops_list[2](ops_list[0](nums[0], ops_list[1](nums[1], nums[2])), nums[3])
            if abs(val2 - 24) < 0.0001:
                return f"({nums[0]}{op_symbols[ops_list[0]]}({nums[1]}{op_symbols[ops_list[1]]}{nums[2]})){op_symbols[ops_list[2]]}{nums[3]}"
            
            # a#(b#(c#d))
            val3 = ops_list[0](nums[0], ops_list[1](nums[1], ops_list[2](nums[2], nums[3])))
            if abs(val3 - 24) < 0.0001:
                return f"{nums[0]}{op_symbols[ops_list[0]]}({nums[1]}{op_symbols[ops_list[1]]}({nums[2]}{op_symbols[ops_list[2]]}{nums[3]})"
            
            # a#((b#c)#d)
            val4 = ops_list[0](nums[0], ops_list[2](ops_list[1](nums[1], nums[2]), nums[3]))
            if abs(val4 - 24) < 0.0001:
                return f"{nums[0]}{op_symbols[ops_list[0]]}(({nums[1]}{op_symbols[ops_list[1]]}{nums[2]}){op_symbols[ops_list[2]]}{nums[3]})"
        except ZeroDivisionError:
            return None
        return None

    for nums in permutations(numbers):
        for op1 in ops:
            for op2 in ops:
                for op3 in ops:
                    result = evaluate(nums, [op1, op2, op3])
                    if result:
                        print(result + " = 24")
                        return

numbers = [1, 1, 4, 6]
generate_expressions(numbers)