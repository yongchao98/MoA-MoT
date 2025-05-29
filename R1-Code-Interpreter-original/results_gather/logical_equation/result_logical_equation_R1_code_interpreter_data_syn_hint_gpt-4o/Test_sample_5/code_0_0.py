# Given numbers
numbers = [1, 3, 4, 9, 16, 27, 36, 80, 121]

# Initialize a dictionary to store the values of each letter
values = {}

# Assign values based on constraints
values['A'] = 9
values['G'] = 36
values['C'] = 36
values['D'] = 9
values['B'] = 27
values['E'] = 121
values['F'] = 80
values['H'] = 3

# The remaining number is for I
values['I'] = list(set(numbers) - set(values.values()))[0]

# Return the values in alphabetical order
result = [values[letter] for letter in sorted(values.keys())]
print(result)