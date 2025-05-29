from itertools import permutations, product

def evaluate_expression(nums, ops, op_symbols):
    # Try different parenthesizations
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"(({nums[0]} {op_symbols[0]} {nums[1]}) {op_symbols[1]} {nums[2]}) {op_symbols[2]} {nums[3]}"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]} {op_symbols[0]} ({nums[1]} {op_symbols[1]} {nums[2]})) {op_symbols[2]} {nums[3]}"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]} {op_symbols[0]} {nums[1]}) {op_symbols[1]} ({nums[2]} {op_symbols[2]} {nums[3]})"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]} {op_symbols[0]} (({nums[1]} {op_symbols[1]} {nums[2]}) {op_symbols[2]} {nums[3]})"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]} {op_symbols[0]} ({nums[1]} {op_symbols[1]} ({nums[2]} {op_symbols[2]} {nums[3]}))"
        
    except ZeroDivisionError:
        pass
    return None

def find_expression_for_24(numbers):
    operations = [lambda x, y: x + y, lambda x, y: x - y, lambda x, y: x * y, lambda x, y: x / y]
    operation_symbols = ['+', '-', '*', '/']
    
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            op_symbols = [operation_symbols[operations.index(op)] for op in ops]
            expression = evaluate_expression(num_perm, ops, op_symbols)
            if expression:
                return f"<<<{expression}>>>"
    return "No solution found."

# Input numbers
numbers = [1, 1, 2, 12]

# Find and print the expression
print(find_expression_for_24(numbers))