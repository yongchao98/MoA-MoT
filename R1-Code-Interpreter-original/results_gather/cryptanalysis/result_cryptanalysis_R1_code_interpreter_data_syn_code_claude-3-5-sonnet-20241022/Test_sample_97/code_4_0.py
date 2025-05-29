def test_password(password):
    def check_guess(guess, feedback):
        g_nums = [int(guess[0]), int(guess[1])]
        g_letters = [guess[2], guess[3]]
        p_nums = [int(password[0]), int(password[1])]
        p_letters = [password[2], password[3]]
        
        if feedback == 1:  # one number correct wrong pos, one small, both letters early
            correct_wrong_pos = sum(1 for i in range(2) for j in range(2) 
                                  if i != j and g_nums[i] == p_nums[j]) == 1
            incorrect_small = sum(1 for i in range(2) 
                                if g_nums[i] < p_nums[i] and g_nums[i] not in p_nums) == 1
            letters_early = all(g_letters[i] < p_letters[i] for i in range(2))
            return correct_wrong_pos and incorrect_small and letters_early
            
        elif feedback == 2:  # both numbers small, one letter correct wrong pos, one early
            nums_small = all(g_nums[i] < p_nums[i] for i in range(2))
            letter_correct_wrong = sum(1 for i in range(2) for j in range(2) 
                                     if i != j and g_letters[i] == p_letters[j]) == 1
            return nums_small and letter_correct_wrong
            
        elif feedback == 3:  # one number correct position, one small, letters early
            correct_pos = sum(1 for i in range(2) if g_nums[i] == p_nums[i]) == 1
            incorrect_small = sum(1 for i in range(2) 
                                if g_nums[i] < p_nums[i] and g_nums[i] != p_nums[i]) == 1
            letters_early = all(g_letters[i] < p_letters[i] for i in range(2))
            return correct_pos and incorrect_small and letters_early
            
        elif feedback == 4:  # one number correct wrong pos, letters incorrect
            correct_wrong_pos = sum(1 for i in range(2) for j in range(2) 
                                  if i != j and g_nums[i] == p_nums[j]) == 1
            return correct_wrong_pos
        
        return True

    # Key test cases that provide the most information
    tests = [
        ("71BQ", 1),  # one correct wrong pos, one small, letters early
        ("40QT", 2),  # both small, one letter correct wrong pos
        ("96MJ", 3),  # one correct pos, one small, letters early
        ("89XJ", 4),  # one correct wrong pos
    ]
    
    # Test each condition
    for guess, feedback in tests:
        if not check_guess(guess, feedback):
            return False
    
    # Additional verification
    p_nums = [int(password[0]), int(password[1])]
    if p_nums[0] < 7 or p_nums[1] < 7:  # Numbers must be large
        return False
    if password[2] < 'Q' or password[3] < 'Q':  # Letters must be late
        return False
    if len(set(password)) != 4:  # No repeating characters
        return False
    
    return True

# Test all possible combinations that meet basic criteria
valid_passwords = []
for n1 in [9]:  # Must be 9 based on "96MJ"
    for n2 in [7, 8]:  # Must be large number
        for l1 in ['Q', 'R', 'S', 'T', 'U', 'V']:
            for l2 in ['Q', 'R', 'S', 'T', 'U', 'V']:
                if l1 != l2:
                    password = [str(n1), str(n2), l1, l2]
                    if test_password(password):
                        valid_passwords.append(password)

print(valid_passwords)