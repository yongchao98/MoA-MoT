def find_24(numbers):
    if len(numbers) == 1:
        return numbers

    results = []
    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                remaining = [numbers[k] for k in range(len(numbers)) if k != i and k != j]
                for op in ['+', '-', '*', '/']:
                    if op == '+':
                        result = numbers[i] + numbers[j]
                        equation = f"({numbers[i]}+{numbers[j]})"
                    elif op == '-':
                        result = numbers[i] - numbers[j]
                        equation = f"({numbers[i]}-{numbers[j]})"
                    elif op == '*':
                        result = numbers[i] * numbers[j]
                        equation = f"({numbers[i]}*{numbers[j]})"
                    elif op == '/':
                        if numbers[j] != 0:
                            result = numbers[i] / numbers[j]
                            equation = f"({numbers[i]}/{numbers[j]})"
                        else:
                            continue

                    new_numbers = remaining + [result]
                    new_equations = find_24(new_numbers)
                    for eq in new_equations:
                        if eq == 24:
                            return [equation]

    return results

numbers = [2, 2, 2, 4]
equations = find_24(numbers)
if equations:
    print(f"<<<{equations[0]}>>>")
else:
    print("No solution found.")