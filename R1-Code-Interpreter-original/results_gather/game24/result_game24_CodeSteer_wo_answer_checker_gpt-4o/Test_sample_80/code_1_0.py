from itertools import permutations

def calculate(a, b, op):
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        return a / b if b != 0 else None

def find_expression(nums, exprs):
    if len(nums) == 1:
        if abs(nums[0] - 24) < 1e-6:  # Check if the result is 24
            return exprs[0]
        else:
            return None

    for i in range(len(nums)):
        for j in range(len(nums)):
            if i != j:
                for op in ['+', '-', '*', '/']:
                    # Create a new list without the i-th and j-th elements
                    new_nums = [nums[k] for k in range(len(nums)) if k != i and k != j]
                    new_exprs = [exprs[k] for k in range(len(exprs)) if k != i and k != j]
                    result = calculate(nums[i], nums[j], op)
                    if result is not None:
                        new_nums.append(result)
                        new_exprs.append(f"({exprs[i]} {op} {exprs[j]})")
                        expr = find_expression(new_nums, new_exprs)
                        if expr:
                            return expr
    return None

numbers = [2, 2, 10, 13]
initial_exprs = [str(num) for num in numbers]

for perm in permutations(numbers):
    expr = find_expression(list(perm), initial_exprs)
    if expr:
        print(f"<<<{expr}>>>")
        break