def check_password(password):
    # Helper function to check if letters are too late in alphabet
    def is_too_late(guess_letter, correct_letter):
        return ord(guess_letter) > ord(correct_letter)
    
    def is_too_early(guess_letter, correct_letter):
        return ord(guess_letter) < ord(correct_letter)

    # Guess 1: 15NT - both numbers wrong, both letters too late
    if '1' in password[:2] or '5' in password[:2]:
        return False
    if 'N' in password[2:] or 'T' in password[2:]:
        return False
    if not (is_too_late('N', password[2]) and is_too_late('T', password[3])):
        return False

    # Guess 2: 21OD - all wrong
    if '2' in password[:2] or '1' in password[:2]:
        return False
    if 'O' in password[2:] or 'D' in password[2:]:
        return False

    # Guess 3: 51BI - all wrong
    if '5' in password[:2] or '1' in password[:2]:
        return False
    if 'B' in password[2:] or 'I' in password[2:]:
        return False

    # Guess 4: 06GE - one number correct in position, one letter correct in position
    correct_count = 0
    if password[0] == '0' or password[1] == '6':
        correct_count += 1
    if password[2] == 'G' or password[3] == 'E':
        correct_count += 1
    if correct_count != 2:
        return False
    if 'E' in password[2:] and is_too_early('E', password[2]):
        return False

    # Guess 5: 32FY - all wrong
    if '3' in password[:2] or '2' in password[:2]:
        return False
    if 'F' in password[2:] or 'Y' in password[2:]:
        return False

    # Guess 6: 13UN - both numbers wrong, both letters too late
    if '1' in password[:2] or '3' in password[:2]:
        return False
    if 'U' in password[2:] or 'N' in password[2:]:
        return False
    if not (is_too_late('U', password[2]) and is_too_late('N', password[3])):
        return False

    # Guess 7: 37DU - all wrong
    if '3' in password[:2] or '7' in password[:2]:
        return False
    if 'D' in password[2:] or 'U' in password[2:]:
        return False

    # Guess 8: 15IG - both numbers wrong, one letter correct but wrong position
    if '1' in password[:2] or '5' in password[:2]:
        return False
    g_correct = ('G' in password[2:] and 'G' != password[2] and 'G' != password[3])
    i_wrong = 'I' not in password[2:]
    if not (g_correct and i_wrong):
        return False

    # Guess 9: 25IA - all wrong
    if '2' in password[:2] or '5' in password[:2]:
        return False
    if 'I' in password[2:] or 'A' in password[2:]:
        return False

    # Guess 10: 93PA - one number correct but wrong position
    if '9' in password[:2] and password[0] != '9' and password[1] != '9':
        pass
    elif '3' in password[:2] and password[0] != '3' and password[1] != '3':
        pass
    else:
        return False
    if 'P' in password[2:] or 'A' in password[2:]:
        return False

    return True

# Generate and test all possible combinations
valid_passwords = []
for n1 in range(10):
    for n2 in range(10):
        if n1 != n2:
            for l1 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                for l2 in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    if l1 != l2:
                        password = f"{n1}{n2}{l1}{l2}"
                        if check_password(password):
                            valid_passwords.append([str(n1), str(n2), l1, l2])

print(valid_passwords)