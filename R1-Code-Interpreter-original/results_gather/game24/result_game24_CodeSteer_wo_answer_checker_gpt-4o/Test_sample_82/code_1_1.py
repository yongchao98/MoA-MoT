from itertools import permutations, product

def evaluate_expression(nums, ops):
    # Evaluate the expression based on the current permutation of numbers and operations
    expression = f"({nums[0]}{ops[0]}{nums[1]}){ops[1]}({nums[2]}{ops[2]}{nums[3]})"
    try:
        if eval(expression) == 24:
            return expression
    except ZeroDivisionError:
        pass
    return None

def find_expression_to_24(numbers):
    # Generate all permutations of the numbers
    for nums in permutations(numbers):
        # Generate all combinations of operations
        for ops in product("+-*/", repeat=3):
            # Evaluate the expression
            result = evaluate_expression(nums, ops)
            if result:
                return result
    return None

numbers = [5, 5, 10, 13]
expression = find_expression_to_24(numbers)

if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")