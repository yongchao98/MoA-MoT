def check_all_conditions(pwd):
    # Must contain B and O (not P or N as they are "too late")
    if 'B' not in pwd or 'O' not in pwd:
        return False
    
    # Guess 2: B must be in different position than in "26BP" (position 2)
    if pwd[2] != 'B':
        return False
    
    # Guess 3:
    # 9 must be in wrong position (not position 1)
    if pwd[1] != '9':
        return False
    # First number must be greater than 3
    if int(pwd[0]) <= 3:
        return False
    # O must be in different position than in "39ON" (position 3)
    if pwd[3] != 'O':
        return False
    
    return True

# Generate possible combinations
first_numbers = ['4', '7', '8']  # Numbers greater than 3
second_number = ['9']  # Must be 9 based on conditions

valid_passwords = []
for n1 in first_numbers:
    pwd = (n1, '9', 'B', 'O')
    if check_all_conditions(pwd):
        valid_passwords.append(pwd)

print(valid_passwords)