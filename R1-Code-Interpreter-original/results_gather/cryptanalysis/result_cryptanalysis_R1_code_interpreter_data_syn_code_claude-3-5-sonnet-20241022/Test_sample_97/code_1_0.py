def check_guess(password, guess, feedback):
    p_nums = [int(password[0]), int(password[1])]
    p_letters = [password[2], password[3]]
    g_nums = [int(guess[0]), int(guess[1])]
    g_letters = [guess[2], guess[3]]
    
    if feedback == "one_correct_wrong_pos_one_small_both_letters_early":
        # Check numbers
        correct_wrong_pos = sum(1 for i in range(2) for j in range(2) if i != j and g_nums[i] == p_nums[j]) == 1
        one_small = sum(1 for i in range(2) if g_nums[i] < p_nums[i]) == 1
        # Check letters
        letters_early = all(g_letters[i] < p_letters[i] for i in range(2))
        return correct_wrong_pos and one_small and letters_early
    
    elif feedback == "both_nums_small_one_letter_correct_wrong_pos_one_early":
        # Check numbers
        nums_small = all(g_nums[i] < p_nums[i] for i in range(2))
        # Check letters
        letter_correct_wrong = sum(1 for i in range(2) for j in range(2) if i != j and g_letters[i] == p_letters[j]) == 1
        return nums_small and letter_correct_wrong
    
    elif feedback == "one_num_correct_pos_one_small":
        correct_pos = sum(1 for i in range(2) if g_nums[i] == p_nums[i]) == 1
        one_small = sum(1 for i in range(2) if g_nums[i] < p_nums[i]) == 1
        return correct_pos and one_small

    return True

def is_valid_password(password):
    guesses = [
        ("71BQ", "one_correct_wrong_pos_one_small_both_letters_early"),
        ("96MJ", "one_num_correct_pos_one_small"),
        ("40QT", "both_nums_small_one_letter_correct_wrong_pos_one_early")
    ]
    
    # Basic rules
    if len(set(password)) != 4:  # All characters must be different
        return False
    
    # Check each guess
    for guess, feedback in guesses:
        if not check_guess(password, guess, feedback):
            return False
    
    return True

# Generate all possible combinations
valid_passwords = []
numbers = range(10)
letters = [chr(i) for i in range(ord('A'), ord('Z')+1)]

for n1 in range(6, 10):  # Numbers must be large based on conditions
    for n2 in range(6, 10):
        if n1 != n2:
            for l1 in letters[16:]:  # Letters must be later in alphabet (Q-Z)
                for l2 in letters[16:]:
                    if l1 != l2:
                        password = [str(n1), str(n2), l1, l2]
                        if is_valid_password(password):
                            valid_passwords.append(password)

print(valid_passwords)