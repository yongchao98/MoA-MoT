def check_all_conditions(pwd):
    # Guess 1: all wrong (already filtered in generation)
    
    # Guess 2: 26BP
    # B or P must be in wrong position (second position)
    if pwd[3] not in ['B', 'P']:
        return False
    
    # Guess 3: 39ON
    # 9 must be in wrong position (second position)
    if pwd[1] != '9':
        return False
    # 3 must be too small compared to first number
    if int(pwd[0]) <= 3:
        return False
    # O or N must be in wrong position (first letter position)
    if pwd[2] not in ['O', 'N']:
        return False
    
    return True

# Generate possible combinations
first_numbers = ['4', '7', '8']  # Numbers greater than 3
second_number = ['9']  # Must be 9 based on conditions
first_letters = ['O', 'N']  # Must be O or N
second_letters = ['B', 'P']  # Must be B or P

valid_passwords = []
for n1 in first_numbers:
    for l1 in first_letters:
        for l2 in second_letters:
            pwd = (n1, '9', l1, l2)
            if check_all_conditions(pwd):
                valid_passwords.append(pwd)

print(valid_passwords)