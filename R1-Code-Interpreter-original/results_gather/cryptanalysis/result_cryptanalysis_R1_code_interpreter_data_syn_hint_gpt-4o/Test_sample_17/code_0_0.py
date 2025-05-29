# Possible numbers and letters based on the deductions
possible_numbers = ['1', '6']
possible_letters = ['Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

# We know 'Q' is in the correct position from Guess 2
# We need to find the correct position for '1' and '6' and the other letter

# Iterate over possible combinations
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:  # Numbers should not repeat
            for letter in possible_letters:
                if letter != 'Q':  # 'Q' is already in the correct position
                    # Check if this combination satisfies all conditions
                    if (num1 == '6' and num2 == '1' and letter > 'L'):
                        print([num1, num2, 'Q', letter])