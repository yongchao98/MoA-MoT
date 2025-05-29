# Possible numbers and letters based on the deductions
possible_numbers = ['6', '1']
possible_letters = ['M', 'O', 'P', 'R', 'S', 'T', 'U', 'V', 'W', 'Z']

# We know 'Q' is in the correct position from Guess 2
# We need to find the correct position for '6' and '1' and the other letter

# Iterate over possible combinations
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:  # Numbers should not repeat
            for letter in possible_letters:
                if letter != 'Q':  # 'Q' is already in the correct position
                    # Check if this combination satisfies all conditions
                    if (num1 == '6' and num2 == '1' and letter > 'L' and letter not in ['N', 'Y', 'X']):
                        print([num1, num2, 'Q', letter])