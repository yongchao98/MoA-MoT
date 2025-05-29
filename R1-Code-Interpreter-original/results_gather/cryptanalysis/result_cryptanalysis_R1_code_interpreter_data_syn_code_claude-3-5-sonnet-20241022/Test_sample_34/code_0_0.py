import string

def is_valid_password(password):
    # Check guess 1: 18VC
    correct_pos = 0
    too_early = 0
    for i in range(2, 4):
        if password[i] == "18VC"[i]:
            correct_pos += 1
        elif password[i] < "18VC"[i]:
            too_early += 1
    if correct_pos != 1 or too_early != 1:
        return False

    # Check guess 2: 54OF
    too_early = 0
    for i in range(2, 4):
        if password[i] == "54OF"[i]:
            return False
        elif password[i] < "54OF"[i]:
            too_early += 1
    if too_early != 2:
        return False

    # Check guess 3: 36AV
    if password[0:2] != "36":
        return False
    correct_wrong_pos = 0
    too_early = 0
    for i in range(2, 4):
        if password[i] == "36AV"[i]:
            return False
        elif "36AV"[i] in password[2:4]:
            correct_wrong_pos += 1
        elif password[i] < "36AV"[i]:
            too_early += 1
    if correct_wrong_pos != 1 or too_early != 1:
        return False

    # Check guess 4: 42GE
    too_early = 0
    for i in range(2, 4):
        if password[i] == "42GE"[i]:
            return False
        elif password[i] < "42GE"[i]:
            too_early += 1
    if too_early != 2:
        return False

    # Check guess 5: 08NQ
    for i in range(2, 4):
        if password[i] == "08NQ"[i]:
            return False

    return True

# Generate all possible combinations
numbers = "36"  # We know these are correct
letters = string.ascii_uppercase
valid_passwords = []

for l1 in letters:
    for l2 in letters:
        if l1 != l2:  # Letters can't repeat
            password = numbers + l1 + l2
            if is_valid_password(password):
                valid_passwords.append([numbers[0], numbers[1], l1, l2])

print(valid_passwords)