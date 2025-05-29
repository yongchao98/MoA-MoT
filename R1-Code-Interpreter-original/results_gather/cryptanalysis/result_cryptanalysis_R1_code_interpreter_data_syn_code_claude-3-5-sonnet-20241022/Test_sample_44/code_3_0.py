def verify_password(password):
    def check_numbers(guess_nums, feedback):
        correct_pos = 0
        correct_val = 0
        for i in range(2):
            if guess_nums[i] == password[i]:
                correct_pos += 1
            elif guess_nums[i] in password[:2]:
                correct_val += 1
        return correct_pos, correct_val

    def check_letters(guess_letters, feedback):
        correct_pos = 0
        correct_val = 0
        for i in range(2):
            if guess_letters[i] == password[i+2]:
                correct_pos += 1
            elif guess_letters[i] in password[2:]:
                correct_val += 1
        return correct_pos, correct_val

    # Test all conditions
    conditions = [
        ("71HT", (1,0), (0,0)),  # one number correct pos, one too small, both letters wrong
        ("19ZS", (0,0), (0,0)),  # all wrong, letters too late
        ("48PA", (0,0), (0,0)),  # all wrong
        ("49FO", (0,0), (0,0)),  # all wrong
        ("25IP", (0,0), (0,0)),  # numbers too small, letters wrong
        ("76KR", (2,0), (0,1)),  # both numbers correct pos, one letter wrong pos
        ("17BY", (0,1), (0,0)),  # one number wrong pos, one too small
        ("31ZP", (0,0), (0,0)),  # numbers too small
        ("02NE", (0,0), (0,0)),  # numbers too small
        ("36ZF", (1,0), (0,0)),  # one number correct pos, one too small
        ("26SR", (1,0), (0,1)),  # one number correct pos, one small, one letter wrong pos
        ("65FG", (0,1), (1,0))   # one number wrong pos, one small, one letter correct pos
    ]

    for guess, num_expected, letter_expected in conditions:
        num_result = check_numbers(guess[:2], None)
        letter_result = check_letters(guess[2:], None)
        
        if num_result != num_expected or letter_result != letter_expected:
            print(f"Failed at {guess}")
            print(f"Number check: Expected {num_expected}, got {num_result}")
            print(f"Letter check: Expected {letter_expected}, got {letter_result}")
            return False
            
    return True

# Test both candidates
passwords = [['7', '6', 'F', 'K'], ['7', '6', 'R', 'G']]
for pwd in passwords:
    print(f"\nTesting {pwd}")
    if verify_password(pwd):
        print("This password satisfies all conditions!")
    else:
        print("This password fails verification.")