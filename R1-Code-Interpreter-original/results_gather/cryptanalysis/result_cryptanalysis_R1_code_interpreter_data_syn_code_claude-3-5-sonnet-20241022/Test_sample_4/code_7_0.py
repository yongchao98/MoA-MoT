def check_guess(pwd):
    pwd_str = ''.join(map(str, pwd))
    
    # Guess 1: 41CK
    def check_first():
        # One number correct in position (must be 4 in first position)
        if pwd[0] != '4':
            return False
        # One number too small (1 must be smaller than second digit)
        if int(pwd[1]) <= 1:
            return False
        # One letter correct but wrong position (K must be present but not in original position)
        if 'K' not in pwd[2:] or pwd[3] == 'K':
            return False
        # One letter incorrect and too early
        if pwd[2] <= 'C' or pwd[3] <= 'C':
            return False
        return True

    # Guess 2: 64DE
    def check_second():
        # One number correct but wrong position (6 must be present but not in original position)
        if '6' not in pwd[:2] or pwd[1] == '6':
            return False
        # Both letters must be after D and E
        if any(l <= 'E' for l in pwd[2:]):
            return False
        return True

    # Guess 3: 87JY
    def check_third():
        # All numbers and letters must be wrong
        if any(x in pwd_str for x in '87JY'):
            return False
        return True

    # Guess 4: 12OD
    def check_fourth():
        # Both numbers must be larger than 1 and 2
        if int(pwd[0]) <= 2 or int(pwd[1]) <= 2:
            return False
        # O must be present but in wrong position
        if 'O' not in pwd[2:] or pwd[3] == 'O':
            return False
        return True

    return all([check_first(), check_second(), check_third(), check_fourth()])

# Test specific combinations based on our deductions:
# 1. First number must be 4
# 2. Second number must be 6
# 3. Must contain K and O
# 4. Letters must be after E
# 5. K and O must be in specific positions to satisfy all conditions

test_cases = [
    [4, 6, 'K', 'O'],
    [4, 6, 'O', 'K'],
    [6, 4, 'K', 'O'],
    [6, 4, 'O', 'K']
]

solutions = []
for test in test_cases:
    if check_guess(test):
        solutions.append(test)

print(solutions)