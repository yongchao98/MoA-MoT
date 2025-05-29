def determine_output(number):
    # Check if the number is even or odd
    if number % 2 == 0:
        return 1
    else:
        return 0

# Test the function with the given examples
examples = [835, 320, 2165, 650]
outputs = [determine_output(num) for num in examples]
print(outputs)

# Determine the output for the input 1715
result = determine_output(1715)
print(result)