# Possible numbers and letters based on the analysis
possible_numbers = [3, 6]
possible_letters = ['Q', 'R']

# Check the conditions
def check_combination(numbers, letters):
    # Condition 6: 96CL
    if numbers[0] == 6 and numbers[1] < 9:
        # Condition 11: 09XR
        if letters[0] == 'Q' and letters[1] == 'R':
            return True
    return False

# Find the correct combination
for num1 in possible_numbers:
    for num2 in possible_numbers:
        if num1 != num2:
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 != letter2:
                        if check_combination([num1, num2], [letter1, letter2]):
                            print([num1, num2, letter1, letter2])