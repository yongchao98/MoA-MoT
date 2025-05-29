from itertools import product

# Define the numbers and target
numbers = [7, 4, 5, 5, 5, 5]
target = 100

# Define possible operations
operations = ['+', '-', '*', '/']

# Function to evaluate an expression safely
def safe_eval(expr):
    try:
        return eval(expr)
    except ZeroDivisionError:
        return None

# Generate all combinations of operations
for ops in product(operations, repeat=5):
    # Generate expression without parentheses
    expr = f"{numbers[0]}{ops[0]}{numbers[1]}{ops[1]}{numbers[2]}{ops[2]}{numbers[3]}{ops[3]}{numbers[4]}{ops[4]}{numbers[5]}"
    if safe_eval(expr) == target:
        print(f"Expression: {expr}")
        break

    # Generate expressions with valid parentheses
    for i in range(5):
        for j in range(i+2, 7):
            if j <= 5:  # Ensure we don't go out of bounds
                expr_list = [str(numbers[0]), ops[0], str(numbers[1]), ops[1], str(numbers[2]), ops[2], str(numbers[3]), ops[3], str(numbers[4]), ops[4], str(numbers[5])]
                expr_list.insert(i * 2, '(')
                expr_list.insert(j * 2 + 1, ')')
                expr = ''.join(expr_list)
                if safe_eval(expr) == target:
                    print(f"Expression: {expr}")
                    break