from itertools import permutations
import operator

def evaluate(nums, ops, parens):
    # Convert string representation to actual operators
    op_dict = {'+': operator.add, '-': operator.sub, '*': operator.mul, '/': operator.truediv}
    
    try:
        if parens == 0:  # ((a#b)#c)#d
            result = op_dict[ops[0]](nums[0], nums[1])
            result = op_dict[ops[1]](result, nums[2])
            result = op_dict[ops[2]](result, nums[3])
        elif parens == 1:  # (a#b)#(c#d)
            result1 = op_dict[ops[0]](nums[0], nums[1])
            result2 = op_dict[ops[2]](nums[2], nums[3])
            result = op_dict[ops[1]](result1, result2)
        elif parens == 2:  # (a#(b#c))#d
            result = op_dict[ops[1]](nums[1], nums[2])
            result = op_dict[ops[0]](nums[0], result)
            result = op_dict[ops[2]](result, nums[3])
        
        return abs(result - 24) < 1e-10
    except:
        return False

def find_24(numbers):
    operators = ['+', '-', '*', '/']
    
    # Try all number permutations
    for nums in permutations(numbers):
        # Try all operator combinations
        for op1 in operators:
            for op2 in operators:
                for op3 in operators:
                    ops = [op1, op2, op3]
                    # Try different parentheses arrangements
                    for parens in range(3):
                        if evaluate(nums, ops, parens):
                            if parens == 0:
                                return f"(({nums[0]}{op1}{nums[1]}){op2}{nums[2]}){op3}{nums[3]}"
                            elif parens == 1:
                                return f"({nums[0]}{op1}{nums[1]}){op2}({nums[2]}{op3}{nums[3]})"
                            else:
                                return f"({nums[0]}{op1}({nums[1]}{op2}{nums[2]})){op3}{nums[3]}"
    return "No solution found"

print(find_24([3, 6, 6, 12]))