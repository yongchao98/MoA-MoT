def check_all_conditions(pwd):
    # Guess 2: 26BP
    # B or P must be in position 3 (wrong position from original guess)
    if pwd[2] not in ['B', 'P']:
        return False
    
    # Guess 3: 39ON
    # 9 must be in position 2 (wrong position)
    if pwd[1] != '9':
        return False
    # First number must be greater than 3
    if int(pwd[0]) <= 3:
        return False
    # O or N must be in position 4 (wrong position from original guess)
    if pwd[3] not in ['O', 'N']:
        return False
    
    return True

# Generate possible combinations
first_numbers = ['4', '7', '8']  # Numbers greater than 3
second_number = ['9']  # Must be 9 based on conditions
first_letters = ['B', 'P']  # Must be B or P in position 3
second_letters = ['O', 'N']  # Must be O or N in position 4

valid_passwords = []
for n1 in first_numbers:
    for l1 in first_letters:
        for l2 in second_letters:
            pwd = (n1, '9', l1, l2)
            if check_all_conditions(pwd):
                valid_passwords.append(pwd)

print(valid_passwords)