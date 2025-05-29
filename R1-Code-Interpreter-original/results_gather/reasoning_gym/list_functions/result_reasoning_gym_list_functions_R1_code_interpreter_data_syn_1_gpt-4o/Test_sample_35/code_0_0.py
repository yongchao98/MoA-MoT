def sum_of_digits(n):
    return sum(int(digit) for digit in str(n))

def determine_output(input_list):
    total_sum = sum(sum_of_digits(num) for num in input_list)
    return [1] if total_sum % 2 != 0 else [0]

input_list = [2, 51, 52, 76, 84, 96]
output = determine_output(input_list)
print(output)