def check_conditions(password):
    n1, n2, l1, l2 = password[0], password[1], password[2], password[3]
    
    # Condition from guess 4 (06GE):
    # - One number (0 or 6) must be correct and in position
    # - One letter (G or E) must be correct and in position
    has_correct_number = (n1 == '0' or n2 == '6')
    has_correct_letter = (l1 == 'G' or l2 == 'E')
    if not (has_correct_number and has_correct_letter):
        return False

    # From guess 8 (15IG):
    # - G must be in the password but in wrong position when guessed at position 3
    if 'G' not in (l1 + l2):
        return False
    
    # From guess 10 (93PA):
    # - Either 9 or 3 must be in the password but in wrong position
    if not ((n1 == '3' and '3' not in guess10_pos) or (n2 == '3' and '3' not in guess10_pos) or
            (n1 == '9' and '9' not in guess10_pos) or (n2 == '9' and '9' not in guess10_pos)):
        return False

    # Numbers that must not appear (from multiple guesses)
    forbidden_numbers = {'1', '2', '5', '7'}
    if n1 in forbidden_numbers or n2 in forbidden_numbers:
        return False

    # Letters that must not appear
    forbidden_letters = {'N', 'T', 'O', 'D', 'B', 'I', 'F', 'Y', 'U', 'P', 'A'}
    if l1 in forbidden_letters or l2 in forbidden_letters:
        return False

    # Letters must be earlier in alphabet than N, T, U
    if not (ord(l1) < ord('N') and ord(l1) < ord('T') and ord(l1) < ord('U') and
            ord(l2) < ord('N') and ord(l2) < ord('T') and ord(l2) < ord('U')):
        return False

    # E must be too early in alphabet when guessed
    if 'E' in (l1 + l2):
        if not (ord(l1) > ord('E') or ord(l2) > ord('E')):
            return False

    return True

# Test all possible combinations
valid_passwords = []
guess10_pos = (0, 1)  # positions where 9 and 3 were guessed in guess 10

for n1 in '0123456789':
    for n2 in '0123456789':
        if n1 != n2:
            for l1 in 'ABCDEFGHIJKLM':  # Only testing letters before N
                for l2 in 'ABCDEFGHIJKLM':
                    if l1 != l2:
                        password = n1 + n2 + l1 + l2
                        if check_conditions(password):
                            valid_passwords.append([n1, n2, l1, l2])

print(valid_passwords)