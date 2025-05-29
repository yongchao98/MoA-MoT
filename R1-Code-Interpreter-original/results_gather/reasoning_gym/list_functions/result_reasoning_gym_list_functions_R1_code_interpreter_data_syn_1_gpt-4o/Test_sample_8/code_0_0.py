def check_divisibility_by_5(number):
    return 1 if number % 5 == 0 else 0

# Test the hypothesis with the given input
input_number = 295
output = check_divisibility_by_5(input_number)
print(output)