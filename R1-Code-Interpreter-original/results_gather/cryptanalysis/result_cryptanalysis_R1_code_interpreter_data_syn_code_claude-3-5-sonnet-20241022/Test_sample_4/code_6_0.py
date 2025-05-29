def check_guess(pwd):
    # Helper function to count matching elements in specific positions
    def count_matches(str1, str2, positions):
        return sum(1 for i in positions if str1[i] == str2[i])
    
    # Helper function to count elements in wrong positions
    def count_wrong_pos(guess, actual, positions):
        return sum(1 for i in positions for j in positions if i != j and guess[i] == actual[j])
    
    # Helper function to count numbers that are too small
    def count_small_numbers(guess, actual, positions):
        return sum(1 for i in positions if int(guess[i]) < int(actual[i]))
    
    # Helper function to count letters that are too early
    def count_early_letters(guess, actual, positions):
        return sum(1 for i in positions if guess[i] < actual[i] and guess[i] not in actual)

    # Guess 1: 41CK
    def check_41CK():
        correct_pos_nums = count_matches('41CK', pwd, [0,1]) == 1
        small_nums = sum(1 for i in [0,1] if int('41'[i]) > int(pwd[i])) == 1
        wrong_pos_letters = count_wrong_pos('41CK', pwd, [2,3]) == 1
        early_letters = sum(1 for i in [2,3] if pwd[i] not in 'CK' and pwd[i] > 'CK'[i-2]) == 1
        return correct_pos_nums and small_nums and wrong_pos_letters and early_letters

    # Guess 2: 64DE
    def check_64DE():
        wrong_pos_nums = count_wrong_pos('64DE', pwd, [0,1]) == 1
        early_letters = sum(1 for i in [2,3] if pwd[i] > 'DE'[i-2]) == 2
        return wrong_pos_nums and early_letters

    # Guess 3: 87JY
    def check_87JY():
        return not any(x in pwd for x in '87JY')

    # Guess 4: 12OD
    def check_12OD():
        small_nums = all(int(pwd[i]) > int('12'[i]) for i in [0,1])
        wrong_pos_letters = count_wrong_pos('12OD', pwd, [2,3]) == 1
        early_letters = sum(1 for i in [2,3] if pwd[i] not in 'OD' and pwd[i] > 'OD'[i-2]) == 1
        return small_nums and wrong_pos_letters and early_letters

    return all([check_41CK(), check_64DE(), check_87JY(), check_12OD()])

# Test all possibilities with known constraints
solutions = []

# We know:
# 1. Numbers must be > 2
# 2. Must contain 4 (from guess 1) and 6 (from guess 2)
# 3. Cannot contain 8,7
# 4. Must contain O (from guess 4)
# 5. Letters must be after D,E

# Try 4 in first position
for l1 in 'KLMNOPQRSTUVWXYZ':
    for l2 in 'KLMNOPQRSTUVWXYZ':
        if l1 != l2 and 'O' in [l1, l2]:
            pwd = f"46{l1}{l2}"
            if check_guess(pwd):
                solutions.append(['4', '6', l1, l2])

# Try 6 in first position
for l1 in 'KLMNOPQRSTUVWXYZ':
    for l2 in 'KLMNOPQRSTUVWXYZ':
        if l1 != l2 and 'O' in [l1, l2]:
            pwd = f"64{l1}{l2}"
            if check_guess(pwd):
                solutions.append(['6', '4', l1, l2])

print(solutions)