# Possible numbers and letters
possible_numbers = [9]
possible_letters = ['X', 'Y', 'Z']

# Check each combination
for number in possible_numbers:
    for letter in possible_letters:
        # Check if the combination satisfies all conditions
        if number > 8 and letter > 'V':
            print([8, number, 'V', letter])