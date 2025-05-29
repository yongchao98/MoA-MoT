input_list = [36, 41, 72, 98, 88, 32, 16, 8, 41]

# Find numbers less than 20
numbers_less_than_20 = [num for num in input_list if num < 20]

# Apply the rule
if numbers_less_than_20:
    result = numbers_less_than_20
else:
    result = [0]

print(result)