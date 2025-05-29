from itertools import permutations, product

numbers = [39, 12, 33, 86]
target = 778
operations = ['+', '-', '*', '/']

def apply_operation(a, b, op):
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        if b != 0 and a % b == 0:  # Ensure division results in an integer
            return a // b
    return None

def find_expression():
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            # Try different parenthesizations
            # (a op1 b) op2 (c op3 d)
            result1 = apply_operation(num_perm[0], num_perm[1], ops[0])
            if result1 is not None:
                result2 = apply_operation(num_perm[2], num_perm[3], ops[2])
                if result2 is not None:
                    final_result = apply_operation(result1, result2, ops[1])
                    if final_result == target:
                        return f"({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]})"
            
            # ((a op1 b) op2 c) op3 d
            result1 = apply_operation(num_perm[0], num_perm[1], ops[0])
            if result1 is not None:
                result2 = apply_operation(result1, num_perm[2], ops[1])
                if result2 is not None:
                    final_result = apply_operation(result2, num_perm[3], ops[2])
                    if final_result == target:
                        return f"(({num_perm[0]} {ops[0]} {num_perm[1]}) {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]}"
            
            # (a op1 (b op2 c)) op3 d
            result1 = apply_operation(num_perm[1], num_perm[2], ops[1])
            if result1 is not None:
                result2 = apply_operation(num_perm[0], result1, ops[0])
                if result2 is not None:
                    final_result = apply_operation(result2, num_perm[3], ops[2])
                    if final_result == target:
                        return f"({num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} {num_perm[2]})) {ops[2]} {num_perm[3]}"
            
            # a op1 ((b op2 c) op3 d)
            result1 = apply_operation(num_perm[1], num_perm[2], ops[1])
            if result1 is not None:
                result2 = apply_operation(result1, num_perm[3], ops[2])
                if result2 is not None:
                    final_result = apply_operation(num_perm[0], result2, ops[0])
                    if final_result == target:
                        return f"{num_perm[0]} {ops[0]} (({num_perm[1]} {ops[1]} {num_perm[2]}) {ops[2]} {num_perm[3]})"
            
            # a op1 (b op2 (c op3 d))
            result1 = apply_operation(num_perm[2], num_perm[3], ops[2])
            if result1 is not None:
                result2 = apply_operation(num_perm[1], result1, ops[1])
                if result2 is not None:
                    final_result = apply_operation(num_perm[0], result2, ops[0])
                    if final_result == target:
                        return f"{num_perm[0]} {ops[0]} ({num_perm[1]} {ops[1]} ({num_perm[2]} {ops[2]} {num_perm[3]}))"
    return "No solution found"

print(find_expression())