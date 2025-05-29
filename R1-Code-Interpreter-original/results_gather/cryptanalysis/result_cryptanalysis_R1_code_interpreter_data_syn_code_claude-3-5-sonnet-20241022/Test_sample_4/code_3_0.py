def check_guess(candidate):
    # Convert candidate to string format if it's a list
    if isinstance(candidate, list):
        candidate = ''.join(map(str, candidate))
    
    # Guess 1: 41CK
    def check_first_guess():
        # One number correct and in position
        correct_pos_num = sum(1 for i in [0,1] if candidate[i] == '41'[i])
        # One number too small
        small_num = sum(1 for i in [0,1] if int('41'[i]) > int(candidate[i]))
        # One letter correct but wrong position
        wrong_pos_letter = sum(1 for i in [2,3] if '41CK'[i] in candidate[2:] and candidate[i] != '41CK'[i])
        # One letter too early
        early_letter = sum(1 for i in [2,3] if '41CK'[i] < candidate[i] and '41CK'[i] not in candidate[2:])
        
        return (correct_pos_num == 1 and small_num == 1 and 
                wrong_pos_letter == 1 and early_letter == 1)

    # Guess 2: 64DE
    def check_second_guess():
        # One number correct but wrong position
        wrong_pos_num = sum(1 for i in [0,1] if '64'[i] in candidate[:2] and candidate[i] != '64'[i])
        # Both letters too early
        early_letters = sum(1 for i in [2,3] if 'DE'[i-2] < candidate[i])
        
        return wrong_pos_num == 1 and early_letters == 2

    # Guess 3: 87JY
    def check_third_guess():
        # Both numbers and letters should be wrong
        return not any(x in candidate for x in '87JY')

    # Guess 4: 12OD
    def check_fourth_guess():
        # Both numbers too small
        small_nums = sum(1 for i in [0,1] if int('12'[i]) < int(candidate[i]))
        # One letter correct but wrong position
        wrong_pos_letter = sum(1 for i in [2,3] if '12OD'[i] in candidate[2:] and candidate[i] != '12OD'[i])
        # One letter too early
        early_letter = sum(1 for i in [2,3] if '12OD'[i] < candidate[i] and '12OD'[i] not in candidate[2:])
        
        return small_nums == 2 and wrong_pos_letter == 1 and early_letter == 1

    return (check_first_guess() and check_second_guess() and 
            check_third_guess() and check_fourth_guess())

# Generate all possible combinations
solutions = []
numbers = range(0, 10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

# Additional constraints:
# - Numbers must be > 2 (from guess 4)
# - Must contain either 4 or 1 in correct position (from guess 1)
# - Must contain 6 (from guess 2)
# - Cannot contain 8,7 (from guess 3)
# - Must contain O (from guess 4)

for n1 in range(3, 10):
    for n2 in range(3, 10):
        if n1 != n2 and not (n1 in [7,8] or n2 in [7,8]):
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:
                        candidate = [str(n1), str(n2), l1, l2]
                        if check_guess(candidate):
                            solutions.append(candidate)

print(solutions)