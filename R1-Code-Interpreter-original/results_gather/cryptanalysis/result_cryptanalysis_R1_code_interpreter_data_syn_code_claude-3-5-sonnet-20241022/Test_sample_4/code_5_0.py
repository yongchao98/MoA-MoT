def check_password(pwd):
    def check_41CK():
        # One number correct and in position
        nums_correct_pos = (pwd[0] == '4' and pwd[1] != '1') or (pwd[0] != '4' and pwd[1] == '1')
        # One number too small
        nums_too_small = (int(pwd[0]) < 4 and int(pwd[1]) >= 1) or (int(pwd[0]) >= 4 and int(pwd[1]) < 1)
        # One letter correct but wrong position
        letters_wrong_pos = ('C' in pwd[2:] and pwd[2] != 'C' and pwd[3] != 'C') or ('K' in pwd[2:] and pwd[2] != 'K' and pwd[3] != 'K')
        # One letter too early
        remaining_letter = 'K' if 'C' in pwd[2:] else 'C' if 'K' in pwd[2:] else None
        if remaining_letter:
            letters_too_early = any(l > remaining_letter for l in pwd[2:] if l != 'C' and l != 'K')
        else:
            letters_too_early = sum(1 for l in pwd[2:] if l > 'C' or l > 'K') == 1
        
        return nums_correct_pos and nums_too_small and letters_wrong_pos and letters_too_early

    def check_64DE():
        # One number correct but wrong position
        nums_wrong_pos = ('6' in pwd[:2] and pwd[0] != '6') or ('4' in pwd[:2] and pwd[1] != '4')
        # Both letters too early in alphabet
        letters_too_early = all(l > 'D' and l > 'E' for l in pwd[2:])
        
        return nums_wrong_pos and letters_too_early

    def check_87JY():
        # Both numbers and letters wrong
        return not any(x in pwd for x in '87JY')

    def check_12OD():
        # Both numbers too small
        nums_too_small = all(int(pwd[i]) > int('12'[i]) for i in range(2))
        # One letter correct but wrong position
        letters_wrong_pos = ('O' in pwd[2:] and pwd[2] != 'O' and pwd[3] != 'O') or ('D' in pwd[2:] and pwd[2] != 'D' and pwd[3] != 'D')
        # One letter too early
        remaining_letter = 'D' if 'O' in pwd[2:] else 'O' if 'D' in pwd[2:] else None
        if remaining_letter:
            letters_too_early = any(l > remaining_letter for l in pwd[2:] if l != 'O' and l != 'D')
        else:
            letters_too_early = sum(1 for l in pwd[2:] if l > 'O' or l > 'D') == 1

        return nums_too_small and letters_wrong_pos and letters_too_early

    return all([check_41CK(), check_64DE(), check_87JY(), check_12OD()])

# Known constraints:
# 1. Numbers must be > 2 (from guess 4)
# 2. Must contain either 4 or 6 (from guesses 1 and 2)
# 3. Cannot contain 8,7 (from guess 3)
# 4. Must contain O (from guess 4)
# 5. Letters must be after D,E (from guess 2)

solutions = []
for n1 in [4, 6]:  # Must contain either 4 or 6
    for n2 in range(3, 10):
        if n2 not in [4, 6, 7, 8] and n1 != n2:  # No repeats, no 7,8
            for l1 in 'KLMNOPQRSTUVWXYZ':  # Letters after E
                for l2 in 'KLMNOPQRSTUVWXYZ':
                    if l1 != l2 and ('O' in [l1, l2]):  # Must contain O
                        pwd = f"{n1}{n2}{l1}{l2}"
                        if check_password(pwd):
                            solutions.append([str(n1), str(n2), l1, l2])

print(solutions)