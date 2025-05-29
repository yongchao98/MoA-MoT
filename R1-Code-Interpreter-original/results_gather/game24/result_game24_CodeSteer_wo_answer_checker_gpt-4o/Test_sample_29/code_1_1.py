from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    a, b, c, d = nums
    op1, op2, op3 = ops
    
    # Try different parenthesizations
    try:
        if eval(f"({a}{op1}{b}){op2}({c}{op3}{d})") == 24:
            return f"(({a}{op1}{b}){op2}({c}{op3}{d}))"
        if eval(f"(({a}{op1}{b}){op2}{c}){op3}{d}") == 24:
            return f"(({a}{op1}{b}){op2}{c}){op3}{d}"
        if eval(f"({a}{op1}({b}{op2}{c})){op3}{d}") == 24:
            return f"({a}{op1}({b}{op2}{c})){op3}{d}"
        if eval(f"{a}{op1}(({b}{op2}{c}){op3}{d})") == 24:
            return f"{a}{op1}(({b}{op2}{c}){op3}{d})"
        if eval(f"{a}{op1}({b}{op2}({c}{op3}{d}))") == 24:
            return f"{a}{op1}({b}{op2}({c}{op3}{d}))"
    except ZeroDivisionError:
        pass
    return None

def find_solution(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product("+-*/", repeat=3):
            result = evaluate_expression(nums, ops)
            if result:
                return f"<<<{result} = 24>>>"
    return "No solution found."

# Input numbers
numbers = [4, 8, 10, 10]

# Find and print the solution
print(find_solution(numbers))