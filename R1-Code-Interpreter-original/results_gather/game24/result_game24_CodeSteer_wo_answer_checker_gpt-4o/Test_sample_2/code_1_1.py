from itertools import permutations

def calculate(a, b, op):
    if op == '+':
        return a + b
    elif op == '-':
        return a - b
    elif op == '*':
        return a * b
    elif op == '/':
        if b != 0:
            return a / b
    return None

def find_24(numbers):
    if len(numbers) == 1:
        if abs(numbers[0] - 24) < 1e-6:  # Check if the result is 24
            return True, str(numbers[0])
        return False, ""

    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                for op in ['+', '-', '*', '/']:
                    # Create a new list without the i-th and j-th elements
                    new_numbers = [numbers[k] for k in range(len(numbers)) if k != i and k != j]
                    result = calculate(numbers[i], numbers[j], op)
                    if result is not None:
                        new_numbers.append(result)
                        found, equation = find_24(new_numbers)
                        if found:
                            return True, f"({numbers[i]} {op} {numbers[j]})" + equation
    return False, ""

# Try all permutations of the numbers
numbers = [1, 1, 2, 8]
for perm in permutations(numbers):
    found, equation = find_24(list(perm))
    if found:
        print(f"<<<{equation} = 24>>>")
        break