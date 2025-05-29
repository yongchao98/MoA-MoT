def apply_constraints():
    # Initial possible values based on feedback
    possible_numbers = set(range(10)) - {0, 3, 4, 5, 9}
    possible_letters = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ') - set('EKOPJC')

    # Further refine based on specific feedback
    # From guess 4: 28JC
    possible_numbers -= {8}  # 8 is too large
    # From guess 5: 98XJ
    possible_numbers -= {9, 8}  # Both numbers are too large
    possible_letters -= {'J'}  # J is too early
    # From guess 6: 37XW
    possible_numbers -= {3, 7}  # Both numbers are incorrect

    # Possible numbers and letters after applying constraints
    possible_numbers = list(possible_numbers)
    possible_letters = list(possible_letters)

    # Debugging output
    print(f"Possible numbers: {possible_numbers}")
    print(f"Possible letters: {possible_letters}")

    return possible_numbers, possible_letters

def find_password():
    possible_numbers, possible_letters = apply_constraints()

    # Try all combinations of the possible numbers and letters
    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 != num2:  # Numbers must not repeat
                for letter1 in possible_letters:
                    for letter2 in possible_letters:
                        if letter1 != letter2:  # Letters must not repeat
                            # Check if the combination satisfies all feedback conditions
                            if (num1 == 2 or num2 == 2) and (letter1 == 'W' or letter2 == 'W') and (letter1 == 'X' or letter2 == 'X'):
                                print(f"Valid combination found: {num1}{num2}{letter1}{letter2}")
                                return [num1, num2, letter1, letter2]
                            else:
                                print(f"Invalid combination: {num1}{num2}{letter1}{letter2}")

# Output the deduced password
password = find_password()
print(f"<<< {password} >>>")