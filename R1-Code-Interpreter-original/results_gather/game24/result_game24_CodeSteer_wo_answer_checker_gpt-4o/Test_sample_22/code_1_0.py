from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Try different groupings of operations
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"((({nums[0]} {ops[0].__name__} {nums[1]}) {ops[1].__name__} {nums[2]}) {ops[2].__name__} {nums[3]}) = 24"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]} {ops[0].__name__} ({nums[1]} {ops[1].__name__} {nums[2]})) {ops[2].__name__} {nums[3]} = 24"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]} {ops[0].__name__} {nums[1]}) {ops[1].__name__} ({nums[2]} {ops[2].__name__} {nums[3]}) = 24"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]} {ops[0].__name__} (({nums[1]} {ops[1].__name__} {nums[2]}) {ops[2].__name__} {nums[3]}) = 24"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]} {ops[0].__name__} ({nums[1]} {ops[1].__name__} ({nums[2]} {ops[2].__name__} {nums[3]})) = 24"
        
    except ZeroDivisionError:
        # Ignore division by zero errors
        pass
    
    return None

def find_solution(numbers):
    # Define the operations
    operations = [lambda x, y: x + y, lambda x, y: x - y, lambda x, y: x * y, lambda x, y: x / y]
    operations_names = ['+', '-', '*', '/']
    
    # Generate all permutations of numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product(operations, repeat=3):
            # Evaluate the expression
            solution = evaluate_expression(nums, ops)
            if solution:
                return f"<<<{solution}>>>"
    
    return "No solution found."

# Input numbers
numbers = [1, 6, 6, 11]

# Find and print the solution
print(find_solution(numbers))