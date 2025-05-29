from itertools import permutations, product

# Define named functions for operations
def add(x, y):
    return x + y

def subtract(x, y):
    return x - y

def multiply(x, y):
    return x * y

def divide(x, y):
    return x / y

# Map operation names for output
operation_names = {add: '+', subtract: '-', multiply: '*', divide: '/'}

def evaluate_expression(nums, ops):
    # Try different parenthesis placements
    try:
        # ((a op1 b) op2 c) op3 d
        result1 = ops[2](ops[1](ops[0](nums[0], nums[1]), nums[2]), nums[3])
        if result1 == 24:
            return f"(({nums[0]}{operation_names[ops[0]]}{nums[1]}){operation_names[ops[1]]}{nums[2]}){operation_names[ops[2]]}{nums[3]}"
        
        # (a op1 (b op2 c)) op3 d
        result2 = ops[2](ops[0](nums[0], ops[1](nums[1], nums[2])), nums[3])
        if result2 == 24:
            return f"({nums[0]}{operation_names[ops[0]]}({nums[1]}{operation_names[ops[1]]}{nums[2]})){operation_names[ops[2]]}{nums[3]}"
        
        # (a op1 b) op2 (c op3 d)
        result3 = ops[1](ops[0](nums[0], nums[1]), ops[2](nums[2], nums[3]))
        if result3 == 24:
            return f"({nums[0]}{operation_names[ops[0]]}{nums[1]}){operation_names[ops[1]]}({nums[2]}{operation_names[ops[2]]}{nums[3]})"
        
        # a op1 ((b op2 c) op3 d)
        result4 = ops[0](nums[0], ops[2](ops[1](nums[1], nums[2]), nums[3]))
        if result4 == 24:
            return f"{nums[0]}{operation_names[ops[0]]}(({nums[1]}{operation_names[ops[1]]}{nums[2]}){operation_names[ops[2]]}{nums[3]})"
        
        # a op1 (b op2 (c op3 d))
        result5 = ops[0](nums[0], ops[1](nums[1], ops[2](nums[2], nums[3])))
        if result5 == 24:
            return f"{nums[0]}{operation_names[ops[0]]}({nums[1]}{operation_names[ops[1]]}({nums[2]}{operation_names[ops[2]]}{nums[3]}))"
        
    except ZeroDivisionError:
        return None

    return None

# Define operations
operations = [add, subtract, multiply, divide]

# Numbers to use
numbers = [1, 3, 7, 8]

# Try all permutations of numbers and operations
for num_perm in permutations(numbers):
    for ops in product(operations, repeat=3):
        expression = evaluate_expression(num_perm, ops)
        if expression:
            print(expression)
            break