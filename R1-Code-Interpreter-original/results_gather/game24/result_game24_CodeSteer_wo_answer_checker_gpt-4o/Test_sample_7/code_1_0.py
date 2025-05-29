from itertools import permutations, product

def apply_operation(a, b, op):
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        if b == 0 or a % b != 0:
            return None
        return a // b

def backtrack(nums, expressions):
    if len(nums) == 1:
        if nums[0] == 24:
            return expressions[0]
        return None

    for i in range(len(nums)):
        for j in range(i + 1, len(nums)):
            for op in '+-*/':
                a, b = nums[i], nums[j]
                result = apply_operation(a, b, op)
                if result is not None:
                    new_nums = nums[:i] + nums[i+1:j] + nums[j+1:] + [result]
                    new_expr = f"({expressions[i]}{op}{expressions[j]})"
                    new_expressions = expressions[:i] + expressions[i+1:j] + expressions[j+1:] + [new_expr]
                    found_expr = backtrack(new_nums, new_expressions)
                    if found_expr:
                        return found_expr
    return None

def find_expression(nums):
    for num_perm in permutations(nums):
        initial_expressions = list(map(str, num_perm))
        result = backtrack(list(num_perm), initial_expressions)
        if result:
            return result
    return None

numbers = [4, 6, 13, 13]
expression = find_expression(numbers)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")