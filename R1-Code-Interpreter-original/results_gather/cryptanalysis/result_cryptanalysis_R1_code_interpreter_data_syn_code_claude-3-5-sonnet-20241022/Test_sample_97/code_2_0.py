def check_condition(password, guess, feedback_type):
    p_nums = [int(password[0]), int(password[1])]
    p_letters = [password[2], password[3]]
    g_nums = [int(guess[0]), int(guess[1])]
    g_letters = [guess[2], guess[3]]
    
    if feedback_type == 1:  # one number correct wrong pos, one small, both letters early
        num_correct_wrong = sum(1 for i in range(2) for j in range(2) if i != j and g_nums[i] == p_nums[j]) == 1
        num_small = sum(1 for i in range(2) if g_nums[i] < p_nums[i]) == 1
        letters_early = all(g_letters[i] < p_letters[i] for i in range(2))
        return num_correct_wrong and num_small and letters_early
        
    elif feedback_type == 2:  # both nums small, one letter correct wrong pos, one early
        nums_small = all(g_nums[i] < p_nums[i] for i in range(2))
        letter_correct_wrong = sum(1 for i in range(2) for j in range(2) if i != j and g_letters[i] == p_letters[j]) == 1
        return nums_small and letter_correct_wrong
        
    elif feedback_type == 3:  # both numbers incorrect, both letters early
        nums_incorrect = all(g_nums[i] != p_nums[j] for i in range(2) for j in range(2))
        letters_early = all(g_letters[i] < p_letters[i] for i in range(2))
        return nums_incorrect and letters_early
        
    elif feedback_type == 4:  # one num correct position, one small, both letters early
        num_correct_pos = sum(1 for i in range(2) if g_nums[i] == p_nums[i]) == 1
        num_small = sum(1 for i in range(2) if g_nums[i] < p_nums[i] and g_nums[i] != p_nums[i]) == 1
        letters_early = all(g_letters[i] < p_letters[i] for i in range(2))
        return num_correct_pos and num_small and letters_early

    return False

def is_valid_password(password):
    conditions = [
        ("71BQ", 1),  # type 1
        ("29BQ", 1),
        ("69BE", 1),
        ("40QT", 2),  # type 2
        ("48NP", 3),  # type 3
        ("15FQ", 3),
        ("83PF", 3),
        ("96MJ", 4),  # type 4
    ]
    
    # Check all conditions
    for guess, feedback_type in conditions:
        if not check_condition(password, guess, feedback_type):
            return False
            
    # Additional checks from other guesses
    if int(password[0]) < 7 or int(password[1]) < 7:  # Numbers must be large
        return False
    if password[2] < 'Q' or password[3] < 'Q':  # Letters must be late in alphabet
        return False
    if len(set(password)) != 4:  # No repeating characters
        return False
    
    return True

# Generate and test all possible combinations
valid_passwords = []
for n1 in range(7, 10):
    for n2 in range(7, 10):
        if n1 != n2:
            for l1 in [chr(i) for i in range(ord('Q'), ord('Z')+1)]:
                for l2 in [chr(i) for i in range(ord('Q'), ord('Z')+1)]:
                    if l1 != l2:
                        password = [str(n1), str(n2), l1, l2]
                        if is_valid_password(password):
                            valid_passwords.append(password)

print(valid_passwords)