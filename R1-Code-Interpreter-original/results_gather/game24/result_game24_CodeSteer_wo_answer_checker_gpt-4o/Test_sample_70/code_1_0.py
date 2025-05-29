def find_24(numbers, equations):
    if len(numbers) == 1:
        if abs(numbers[0] - 24) < 1e-6:  # Check if the result is 24
            return equations[0]
        else:
            return None

    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]
                for op in ['+', '-', '*', '/']:
                    if op == '+':
                        result = numbers[i] + numbers[j]
                        equation = f"({equations[i]}+{equations[j]})"
                    elif op == '-':
                        result = numbers[i] - numbers[j]
                        equation = f"({equations[i]}-{equations[j]})"
                    elif op == '*':
                        result = numbers[i] * numbers[j]
                        equation = f"({equations[i]}*{equations[j]})"
                    elif op == '/':
                        if numbers[j] != 0:
                            result = numbers[i] / numbers[j]
                            equation = f"({equations[i]}/{equations[j]})"
                        else:
                            continue

                    new_numbers = remaining + [result]
                    new_equations = [equations[k] for k in range(len(numbers)) if k != i and k != j] + [equation]
                    solution = find_24(new_numbers, new_equations)
                    if solution:
                        return solution

    return None

numbers = [2, 2, 2, 4]
equations = [str(num) for num in numbers]
solution = find_24(numbers, equations)
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")