def validate_combination(n1, n2, l1, l2):
    # Convert to strings for easier comparison
    pwd = [str(n1), str(n2), l1, l2]
    
    # Guess 1: 41CK
    def check_guess1():
        # Count correct numbers in position
        correct_num_pos = (pwd[0] == '4')  # Must be 4 in first position
        # Count numbers too small
        small_num = (int(pwd[1]) > 1)  # Second number must be larger than 1
        # Count correct letters in wrong position
        k_present = ('K' in [l1, l2])  # K must be present
        # Count early letters
        letters_after_c = (l1 > 'C' and l2 > 'C')  # Both letters must be after C
        return correct_num_pos and small_num and k_present and letters_after_c

    # Guess 2: 64DE
    def check_guess2():
        # One number correct but wrong position
        six_present = ('6' in [pwd[0], pwd[1]])  # 6 must be present
        # Both letters too early
        letters_after_e = (l1 > 'E' and l2 > 'E')  # Both letters must be after E
        return six_present and letters_after_e

    # Guess 3: 87JY
    def check_guess3():
        # No 8, 7, J, or Y should be present
        return not any(x in pwd for x in ['8', '7', 'J', 'Y'])

    # Guess 4: 12OD
    def check_guess4():
        # Both numbers too small
        numbers_larger = (int(pwd[0]) > 2 and int(pwd[1]) > 2)
        # O must be present
        o_present = ('O' in [l1, l2])
        return numbers_larger and o_present

    return all([
        check_guess1(),
        check_guess2(),
        check_guess3(),
        check_guess4()
    ])

# We know:
# 1. First number must be 4
# 2. Second number must be 6
# 3. Must contain both K and O
# 4. Letters must be after E
# 5. Cannot contain 8,7,J,Y

# Test only the most likely combinations based on our deductions
test_cases = [
    (4, 6, 'K', 'O'),
    (4, 6, 'O', 'K'),
    (4, 6, 'M', 'O'),
    (4, 6, 'O', 'M'),
    (4, 6, 'K', 'P'),
    (4, 6, 'P', 'K')
]

solutions = []
for n1, n2, l1, l2 in test_cases:
    if validate_combination(n1, n2, l1, l2):
        solutions.append([str(n1), str(n2), l1, l2])

print(solutions)