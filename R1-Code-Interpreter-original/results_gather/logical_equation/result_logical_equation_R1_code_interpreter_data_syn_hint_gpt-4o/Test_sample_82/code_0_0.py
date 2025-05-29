# Assign the values to the letters
values = {
    'A': 50,
    'B': None,  # B is not involved in any equation, so it must be the remaining number
    'C': 24,
    'D': 2,
    'E': 45,
    'F': 10,
    'G': 3,
    'H': 7,
    'I': 28,
    'J': 96,
    'K': 16,
    'L': 15,
    'M': 5
}

# List of possible values
possible_values = [2, 3, 5, 7, 10, 15, 16, 24, 28, 36, 45, 50, 96]

# Find the value for B
for value in possible_values:
    if value not in values.values():
        values['B'] = value
        break

# Return the values in alphabetical order
result = [values[letter] for letter in sorted(values.keys())]
print(result)