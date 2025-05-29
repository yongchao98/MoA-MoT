from itertools import permutations, product

def evaluate_expression(nums, ops, op_symbols):
    # Try different parenthesis placements
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"(({nums[0]} {op_symbols[ops[0]]} {nums[1]}) {op_symbols[ops[1]]} {nums[2]}) {op_symbols[ops[2]]} {nums[3]}"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]} {op_symbols[ops[0]]} ({nums[1]} {op_symbols[ops[1]]} {nums[2]})) {op_symbols[ops[2]]} {nums[3]}"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]} {op_symbols[ops[0]]} {nums[1]}) {op_symbols[ops[1]]} ({nums[2]} {op_symbols[ops[2]]} {nums[3]})"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]} {op_symbols[ops[0]]} (({nums[1]} {op_symbols[ops[1]]} {nums[2]}) {op_symbols[ops[2]]} {nums[3]})"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]} {op_symbols[ops[0]]} ({nums[1]} {op_symbols[ops[1]]} ({nums[2]} {op_symbols[ops[2]]} {nums[3]}))"
        
    except ZeroDivisionError:
        pass
    return None

def find_expression_for_24(numbers):
    operations = [
        (lambda x, y: x + y, '+'),
        (lambda x, y: x - y, '-'),
        (lambda x, y: x * y, '*'),
        (lambda x, y: x / y, '/')
    ]
    op_dict = {op[0]: op[1] for op in operations}
    for num_perm in permutations(numbers):
        for ops in product([op[0] for op in operations], repeat=3):
            expression = evaluate_expression(num_perm, ops, op_dict)
            if expression:
                return expression
    return "No solution found"

numbers = [3, 5, 5, 6]
solution = find_expression_for_24(numbers)
print(solution)