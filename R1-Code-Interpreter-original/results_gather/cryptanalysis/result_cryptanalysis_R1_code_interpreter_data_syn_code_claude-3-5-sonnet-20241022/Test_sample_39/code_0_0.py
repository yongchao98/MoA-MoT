from itertools import permutations

def check_guess(password, guess, feedback):
    # Unpack feedback
    if guess == "87ID":
        num_correct_wrong_pos = 1
        num_too_large = 1
        letter_correct_wrong_pos = 1
        letter_too_early = 1
    else:  # 50FH
        num_correct_wrong_pos = 0
        num_too_large = 0
        letter_correct_wrong_pos = 1
        letter_too_early = 1
    
    # Check numbers
    correct_num_wrong_pos = 0
    too_large_nums = 0
    too_small_nums = 0
    for i in range(2):
        if guess[i] in password[:2] and guess[i] != password[i]:
            correct_num_wrong_pos += 1
        elif int(guess[i]) > int(password[i]):
            too_large_nums += 1
        elif int(guess[i]) < int(password[i]):
            too_small_nums += 1
    
    # Check letters
    correct_letter_wrong_pos = 0
    too_early_letters = 0
    for i in range(2,4):
        if guess[i] in password[2:] and guess[i] != password[i]:
            correct_letter_wrong_pos += 1
        elif guess[i] < password[i]:
            too_early_letters += 1
    
    # Check all conditions
    if guess == "87ID":
        return (correct_num_wrong_pos == num_correct_wrong_pos and
                too_large_nums == num_too_large and
                correct_letter_wrong_pos == letter_correct_wrong_pos and
                too_early_letters == letter_too_early)
    else:  # 50FH
        return (too_small_nums == 2 and
                correct_letter_wrong_pos == letter_correct_wrong_pos and
                too_early_letters == letter_too_early)

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

valid_passwords = []

# Try all possible combinations
for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # no repeating numbers
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # no repeating letters
                        password = n1 + n2 + l1 + l2
                        if (check_guess(password, "87ID", None) and 
                            check_guess(password, "50FH", None)):
                            valid_passwords.append([n1, n2, l1, l2])

print(valid_passwords)