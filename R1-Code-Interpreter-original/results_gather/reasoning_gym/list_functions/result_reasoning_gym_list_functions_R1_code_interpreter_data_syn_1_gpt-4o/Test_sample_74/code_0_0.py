def determine_output(number):
    # Check divisibility by 5
    if number % 5 == 0:
        return 0
    else:
        return 1

# Test the function with the given input
input_number = 1475
output = determine_output(input_number)
print(output)