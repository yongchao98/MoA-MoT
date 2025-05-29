from itertools import permutations, product

numbers = [3, 8, 8, 9]
target = 24
operations = ['+', '-', '*', '/']

def evaluate_expression(nums, ops):
    try:
        result = nums[0]
        for i in range(3):
            if ops[i] == '+':
                result += nums[i + 1]
            elif ops[i] == '-':
                result -= nums[i + 1]
            elif ops[i] == '*':
                result *= nums[i + 1]
            elif ops[i] == '/':
                result /= nums[i + 1]
        return result
    except ZeroDivisionError:
        return None

def find_solution(numbers, target):
    for num_perm in permutations(numbers):
        for ops in product(operations, repeat=3):
            if evaluate_expression(num_perm, ops) == target:
                return f"({num_perm[0]}{ops[0]}{num_perm[1]}){ops[1]}({num_perm[2]}{ops[2]}{num_perm[3]})"
    return "No solution found"

solution = find_solution(numbers, target)
print(solution)