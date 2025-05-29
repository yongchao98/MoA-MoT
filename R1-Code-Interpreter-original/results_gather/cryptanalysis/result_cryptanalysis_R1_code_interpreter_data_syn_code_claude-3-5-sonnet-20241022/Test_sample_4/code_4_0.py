def validate_password(password):
    # Helper function to check if a letter is too early in alphabet
    def is_letter_early(guess_letter, actual_letter):
        return guess_letter < actual_letter
    
    # Guess 1: 41CK
    def check_guess1():
        # One number correct and in position
        correct_pos_num = sum(1 for i in range(2) if str(password[i]) == '41'[i])
        # One number too small
        small_num = sum(1 for i in range(2) if int('41'[i]) > int(password[i]))
        # One letter correct but wrong position
        wrong_pos_letter = sum(1 for i in [2,3] for j in [2,3] if i != j and password[i] == '41CK'[j])
        # One letter too early
        early_letter = sum(1 for i in [2,3] if password[i] not in 'CK' and is_letter_early('CK'[i-2], password[i]))
        
        return (correct_pos_num == 1 and small_num == 1 and 
                wrong_pos_letter == 1 and early_letter == 1)

    # Guess 2: 64DE
    def check_guess2():
        # One number correct but wrong position
        wrong_pos_num = sum(1 for i in range(2) for j in range(2) if i != j and str(password[i]) == '64'[j])
        # Both letters too early
        early_letters = sum(1 for i in [2,3] if is_letter_early('DE'[i-2], password[i]))
        
        return wrong_pos_num == 1 and early_letters == 2

    # Guess 3: 87JY
    def check_guess3():
        # All numbers and letters should be wrong
        nums_wrong = all(str(n) not in '87' for n in password[:2])
        letters_wrong = all(l not in 'JY' for l in password[2:])
        return nums_wrong and letters_wrong

    # Guess 4: 12OD
    def check_guess4():
        # Both numbers too small
        small_nums = sum(1 for i in range(2) if int('12'[i]) < int(password[i]))
        # One letter correct but wrong position
        wrong_pos_letter = sum(1 for i in [2,3] for j in [2,3] if i != j and password[i] == '12OD'[j])
        # One letter too early
        early_letter = sum(1 for i in [2,3] if password[i] not in 'OD' and is_letter_early('OD'[i-2], password[i]))
        
        return small_nums == 2 and wrong_pos_letter == 1 and early_letter == 1

    return (check_guess1() and check_guess2() and check_guess3() and check_guess4())

# Generate and test combinations with known constraints
solutions = []

# We know:
# - Numbers must be > 2 (from guess 4)
# - Must contain 4 (likely in first position from guess 1)
# - Must contain 6 (from guess 2)
# - Cannot contain 8,7 (from guess 3)
# - Must contain O (from guess 4)

for n2 in range(3, 10):
    if n2 != 4 and n2 not in [7, 8]:  # Exclude 4 (used as n1), 7, and 8
        for l1 in 'KLMNPQRSTUVWXYZ':  # Letters after D (from guess 2)
            for l2 in 'KLMNPQRSTUVWXYZ':
                if l1 != l2 and 'O' in [l1, l2]:  # Must contain O
                    test = [4, n2, l1, l2]  # 4 in first position
                    if validate_password(test):
                        solutions.append(test)

# Try with 6 in first position as well
for n2 in range(3, 10):
    if n2 != 6 and n2 not in [7, 8]:  # Exclude 6 (used as n1), 7, and 8
        for l1 in 'KLMNPQRSTUVWXYZ':  # Letters after D
            for l2 in 'KLMNPQRSTUVWXYZ':
                if l1 != l2 and 'O' in [l1, l2]:  # Must contain O
                    test = [6, n2, l1, l2]  # 6 in first position
                    if validate_password(test):
                        solutions.append(test)

print(solutions)