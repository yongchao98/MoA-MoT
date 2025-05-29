def find_number_with_highest_tens_digit(numbers):
    # Sort numbers by tens digit and then by value
    sorted_numbers = sorted(numbers, key=lambda x: (x // 10, x), reverse=True)
    # Return the first number in the sorted list
    return sorted_numbers[0]

input_numbers = [37, 26, 58, 11, 43, 16]
result = find_number_with_highest_tens_digit(input_numbers)
print(result)