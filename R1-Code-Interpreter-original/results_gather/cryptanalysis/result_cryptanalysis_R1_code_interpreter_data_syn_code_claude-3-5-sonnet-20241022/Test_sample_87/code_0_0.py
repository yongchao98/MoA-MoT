def check_guess(guess, actual):
    # Check if a guess matches all the feedback conditions
    def check_numbers(g_nums, a_nums):
        correct_pos = sum(1 for i in range(2) if g_nums[i] == a_nums[i])
        correct_wrong_pos = sum(1 for i in range(2) if g_nums[i] in a_nums and g_nums[i] != a_nums[i])
        return correct_pos, correct_wrong_pos

    def check_letters(g_lets, a_lets):
        correct_pos = sum(1 for i in range(2) if g_lets[i] == a_lets[i])
        correct_wrong_pos = sum(1 for i in range(2) if g_lets[i] in a_lets and g_lets[i] != a_lets[i])
        return correct_pos, correct_wrong_pos

    g_nums = [int(guess[0]), int(guess[1])]
    g_lets = [guess[2], guess[3]]
    a_nums = [int(actual[0]), int(actual[1])]
    a_lets = [actual[2], actual[3]]
    
    return check_numbers(g_nums, a_nums), check_letters(g_lets, a_lets)

def is_valid_password(password):
    # Convert password to proper format for checking
    password_str = f"{password[0]}{password[1]}{password[2]}{password[3]}"
    
    # Check all conditions
    guesses = [
        ("27JR", ((0,0), (0,0))),  # both numbers and letters incorrect
        ("67GP", ((0,0), (0,1))),  # both numbers too large, one letter wrong pos
        ("09CO", ((0,1), (0,0))),  # one number wrong pos, one too large
        ("57VJ", ((0,0), (0,0))),  # both numbers too large
        ("06HO", ((0,1), (0,0))),  # one number wrong pos, one too large
        ("67JQ", ((0,0), (0,1))),  # both numbers too large, one letter wrong pos
        ("34TS", ((1,0), (0,0)))   # one number correct pos, one too large
    ]
    
    for guess, expected_result in guesses:
        if check_guess(guess, password_str) != expected_result:
            return False
            
    return True

# Generate all possible combinations
valid_passwords = []
for n1 in range(5):  # numbers must be less than 5
    for n2 in range(5):
        if n1 == n2:
            continue
        for l1 in 'ABCDEFGHIJKLMNOPQR':  # letters before S
            for l2 in 'ABCDEFGHIJKLMNOPQR':
                if l1 == l2:
                    continue
                password = [str(n1), str(n2), l1, l2]
                if is_valid_password(password):
                    valid_passwords.append(password)

print(valid_passwords)