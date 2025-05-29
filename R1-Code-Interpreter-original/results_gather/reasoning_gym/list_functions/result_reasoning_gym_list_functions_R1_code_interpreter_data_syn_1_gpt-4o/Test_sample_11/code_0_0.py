def determine_output(number):
    if number % 5 == 0:
        return 1
    else:
        return 0

input_number = 2265
output = determine_output(input_number)
print(output)