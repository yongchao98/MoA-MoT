def check_conditions(password):
    # Convert password to list of strings
    if isinstance(password, list):
        pass_list = password
    else:
        pass_list = list(str(password))
    
    # Condition 1: 09OW - both numbers correct position, both letters wrong
    if pass_list[0] != '0' or pass_list[1] != '9':
        return False
    if pass_list[2] == 'O' or pass_list[3] == 'W':
        return False
    
    # Condition 2: 71BK - all wrong
    if '7' in pass_list or '1' in pass_list or 'B' in pass_list or 'K' in pass_list:
        return False
    
    # Condition 3: 40TV - one number wrong position, T&V too late
    count_40 = sum(1 for x in ['4', '0'] if x in pass_list)
    if count_40 != 1:
        return False
    if 'T' in pass_list or 'V' in pass_list:
        return False
    
    # Condition 4: 92SU - one number wrong position, S correct position, U too late
    count_92 = sum(1 for x in ['9', '2'] if x in pass_list)
    if count_92 != 1:
        return False
    if pass_list[2] != 'S':
        return False
    if 'U' in pass_list:
        return False
    
    # Condition 5: 41SD - numbers wrong, S&D correct positions
    if '4' in pass_list or '1' in pass_list:
        return False
    if pass_list[2] != 'S' or pass_list[3] != 'D':
        return False
    
    return True

# Generate all possible combinations
import string
valid_passwords = []

# We know first two digits are 0 and 9
numbers = ['0', '9']
letters = list(string.ascii_uppercase)

# Remove letters we know are wrong or too late
letters = [l for l in letters if l not in ['O', 'W', 'B', 'K', 'T', 'V', 'U']]

# We know S and D are in positions 2 and 3
if check_conditions(['0', '9', 'S', 'D']):
    valid_passwords.append(['0', '9', 'S', 'D'])

print(valid_passwords)