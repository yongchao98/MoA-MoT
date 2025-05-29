from itertools import permutations
import operator

# Define the operations
ops = [operator.add, operator.sub, operator.mul, operator.truediv]
ops_symbols = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    # Try all possible ways to insert parentheses
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"((({nums[0]} {ops_symbols[0]} {nums[1]}) {ops_symbols[1]} {nums[2]}) {ops_symbols[2]} {nums[3]})"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]} {ops_symbols[0]} ({nums[1]} {ops_symbols[1]} {nums[2]})) {ops_symbols[2]} {nums[3]}"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]} {ops_symbols[0]} {nums[1]}) {ops_symbols[1]} ({nums[2]} {ops_symbols[2]} {nums[3]})"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]} {ops_symbols[0]} (({nums[1]} {ops_symbols[1]} {nums[2]}) {ops_symbols[2]} {nums[3]})"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]} {ops_symbols[0]} ({nums[1]} {ops_symbols[1]} ({nums[2]} {ops_symbols[2]} {nums[3]}))"
        
    except ZeroDivisionError:
        pass
    return None

def find_expression(nums):
    for num_perm in permutations(nums):
        for op_perm in permutations(ops, 3):
            expression = evaluate_expression(num_perm, op_perm)
            if expression:
                return expression
    return "No solution found"

# Input numbers
numbers = [1, 1, 5, 5]

# Find and print the expression
expression = find_expression(numbers)
print(f"<<<{expression}>>>")