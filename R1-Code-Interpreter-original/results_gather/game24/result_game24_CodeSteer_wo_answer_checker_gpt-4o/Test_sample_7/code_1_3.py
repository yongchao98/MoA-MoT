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

def backtrack(nums, ops, expression):
    if len(nums) == 1:
        if nums[0] == 24:
            return expression
        return None

    for i in range(len(nums)):
        for j in range(i + 1, len(nums)):
            for op in ops:
                a, b = nums[i], nums[j]
                result = apply_operation(a, b, op)
                if result is not None:
                    new_nums = nums[:i] + nums[i+1:j] + nums[j+1:] + [result]
                    new_expr = f"({expression[i]}{op}{expression[j]})"
                    new_expression = expression[:i] + expression[i+1:j] + expression[j+1:] + [new_expr]
                    found_expr = backtrack(new_nums, ops, new_expression)
                    if found_expr:
                        return found_expr
    return None

def find_expression(nums):
    ops = '+-*/'
    for num_perm in permutations(nums):
        initial_expr = list(map(str, num_perm))
        result = backtrack(list(num_perm), ops, initial_expr)
        if result:
            return result
    return None

numbers = [4, 6, 13, 13]
expression = find_expression(numbers)
if expression:
    print(f"<<<{expression}>>>")
else:
    print("No solution found.")