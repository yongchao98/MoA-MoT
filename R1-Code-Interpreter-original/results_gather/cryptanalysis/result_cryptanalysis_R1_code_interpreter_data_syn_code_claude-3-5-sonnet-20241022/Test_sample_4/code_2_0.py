def check_conditions(test):
    # Convert test to string if it's a list
    if isinstance(test, list):
        test = ''.join(map(str, test))
    
    # Check against 41CK
    def check_41CK():
        # We know 4 must be in first position
        if test[0] != '4':
            return False
        small_nums = sum(1 for i in [0,1] if int('41'[i]) < int(test[i]))
        correct_wrong_pos_letters = sum(1 for i in [2,3] if '41CK'[i] in test[2:] and '41CK'[i] != test[i])
        early_letters = sum(1 for i in [2,3] if '41CK'[i] < test[i] and '41CK'[i] not in test[2:])
        return small_nums == 1 and correct_wrong_pos_letters == 1 and early_letters == 1

    # Check against 64DE
    def check_64DE():
        # Must contain 6 in wrong position
        if '6' not in test[:2] or test.index('6') == 1:  # 6 must be present but not in second position
            return False
        early_letters = sum(1 for i in [2,3] if 'DE'[i-2] < test[i])
        return early_letters == 2

    # Check against 87JY
    def check_87JY():
        nums_wrong = all(n not in test[:2] for n in '87')
        letters_wrong = all(l not in test[2:] for l in 'JY')
        return nums_wrong and letters_wrong

    # Check against 12OD
    def check_12OD():
        small_nums = sum(1 for i in [0,1] if int('12'[i]) < int(test[i]))
        # O must be present but in wrong position
        if 'O' not in test[2:] or test.index('O') == 3:
            return False
        early_letters = sum(1 for i in [2,3] if '12OD'[i] < test[i] and '12OD'[i] not in test[2:])
        return small_nums == 2 and early_letters == 1

    return (check_41CK() and check_64DE() and check_87JY() and check_12OD())

# Generate and test all possible combinations
solutions = []
for n2 in range(3, 10):  # numbers must be > 2
    if n2 != 4:  # first number must be 4
        test = ['4', str(n2), 'K', 'O']
        if '6' in [test[0], test[1]]:  # must contain 6
            if check_conditions(test):
                solutions.append(test)

print(solutions)