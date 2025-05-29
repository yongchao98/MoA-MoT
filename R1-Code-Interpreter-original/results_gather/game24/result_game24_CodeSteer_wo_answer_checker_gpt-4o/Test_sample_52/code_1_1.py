import itertools
import operator

def evaluate_expression(nums, ops):
    # Try different parenthesizations
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"(({nums[0]}{ops[0].__name__}{nums[1]}){ops[1].__name__}{nums[2]}){ops[2].__name__}{nums[3]}"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]}{ops[0].__name__}({nums[1]}{ops[1].__name__}{nums[2]})){ops[2].__name__}{nums[3]}"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]}{ops[0].__name__}{nums[1]}){ops[1].__name__}({nums[2]}{ops[2].__name__}{nums[3]})"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]}{ops[0].__name__}(({nums[1]}{ops[1].__name__}{nums[2]}){ops[2].__name__}{nums[3]})"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]}{ops[0].__name__}({nums[1]}{ops[1].__name__}({nums[2]}{ops[2].__name__}{nums[3]}))"
        
    except ZeroDivisionError:
        pass
    return None

def find_expression_to_24(numbers):
    operations = [operator.add, operator.sub, operator.mul, operator.truediv]
    for nums in itertools.permutations(numbers):
        for ops in itertools.product(operations, repeat=3):
            expression = evaluate_expression(nums, ops)
            if expression:
                return f"<<<{expression}>>>"
    return "No solution found."

numbers = [4, 5, 9, 9]
print(find_expression_to_24(numbers))