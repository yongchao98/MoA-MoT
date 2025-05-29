def find_24(numbers, equation=""):
    if len(numbers) == 1:
        if abs(numbers[0] - 24) < 1e-6:  # Check if the number is close to 24
            return equation
        else:
            return None

    for i in range(len(numbers)):
        for j in range(len(numbers)):
            if i != j:
                # Create a new list without the i-th and j-th elements
                new_numbers = [numbers[k] for k in range(len(numbers)) if k != i and k != j]

                # Try all operations
                for op in ['+', '-', '*', '/']:
                    if op == '+':
                        result = numbers[i] + numbers[j]
                        new_equation = f"({equation}{numbers[i]}+{numbers[j]})"
                    elif op == '-':
                        result = numbers[i] - numbers[j]
                        new_equation = f"({equation}{numbers[i]}-{numbers[j]})"
                    elif op == '*':
                        result = numbers[i] * numbers[j]
                        new_equation = f"({equation}{numbers[i]}*{numbers[j]})"
                    elif op == '/':
                        if numbers[j] != 0:
                            result = numbers[i] / numbers[j]
                            new_equation = f"({equation}{numbers[i]}/{numbers[j]})"
                        else:
                            continue

                    # Add the result to the new list and recurse
                    new_numbers_with_result = new_numbers + [result]
                    sub_equation = find_24(new_numbers_with_result, new_equation)
                    if sub_equation:
                        return sub_equation

    return None

numbers = [7, 7, 11, 12]
solution = find_24(numbers)
if solution:
    print(f"<<<{solution}>>>")
else:
    print("No solution found.")