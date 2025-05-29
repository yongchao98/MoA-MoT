import itertools
import operator

# Define the operations
ops = {
    '+': operator.add,
    '-': operator.sub,
    '*': operator.mul,
    '/': operator.truediv
}

# Function to evaluate an expression with given numbers and operations
def evaluate_expression(numbers, operations):
    try:
        # Try different parenthesizations
        for a, b, c in itertools.permutations(range(3), 3):
            # ((a op b) op c)
            if ops[operations[0]](ops[operations[1]](numbers[a], numbers[b]), numbers[c]) == 24:
                return f"(({numbers[a]} {operations[1]} {numbers[b]}) {operations[0]} {numbers[c]})"
            # (a op (b op c))
            if ops[operations[0]](numbers[a], ops[operations[1]](numbers[b], numbers[c])) == 24:
                return f"({numbers[a]} {operations[0]} ({numbers[b]} {operations[1]} {numbers[c]}))"
    except ZeroDivisionError:
        pass
    return None

# Main function to find the expression
def find_expression_to_24(numbers):
    # Generate all permutations of numbers
    for num_perm in itertools.permutations(numbers):
        # Generate all combinations of operations
        for ops_comb in itertools.product(ops.keys(), repeat=2):
            # Evaluate the expression
            result = evaluate_expression(num_perm, ops_comb)
            if result:
                return f"<<<{result}>>>"
    return "No solution found."

# Example input
numbers = [1, 1, 2, 13]

# Find and print the expression
print(find_expression_to_24(numbers))