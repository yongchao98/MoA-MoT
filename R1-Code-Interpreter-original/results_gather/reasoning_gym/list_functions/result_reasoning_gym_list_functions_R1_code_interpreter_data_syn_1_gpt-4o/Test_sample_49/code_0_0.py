def check_property(numbers):
    # Check if there is any even number in the list
    for number in numbers:
        if number % 2 == 0:
            return 1
    return 0

# Test the function with the given input
input_list = [3, 30, 32, 58, 58, 72, 84, 90, 99]
output = check_property(input_list)
print([output])