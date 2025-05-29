def verify_guess(password, guess, feedback):
    # Split password and guess into numbers and letters
    pass_nums = password[:2]
    pass_letters = password[2:]
    guess_nums = guess[:2]
    guess_letters = guess[2:]
    
    # Check numbers
    correct_num_pos = sum(1 for i in range(2) if guess_nums[i] == pass_nums[i])
    correct_num_wrong_pos = sum(1 for x in guess_nums if x in pass_nums) - correct_num_pos
    
    # Check letters
    correct_letter_pos = sum(1 for i in range(2) if guess_letters[i] == pass_letters[i])
    correct_letter_wrong_pos = sum(1 for x in guess_letters if x in pass_letters) - correct_letter_pos
    
    return correct_num_pos, correct_num_wrong_pos, correct_letter_pos, correct_letter_wrong_pos

def test_password(password):
    # Test all conditions
    guesses_and_feedback = [
        ("71HT", (1,0,0,0)),  # one number correct position, one too small
        ("19ZS", (0,0,0,0)),  # both numbers incorrect, letters too late
        ("48PA", (0,0,0,0)),  # all incorrect
        ("49FO", (0,0,0,0)),  # all incorrect
        ("25IP", (0,0,0,0)),  # numbers too small
        ("76KR", (2,0,0,1)),  # both numbers correct, one letter wrong pos
        ("17BY", (0,1,0,0)),  # one number wrong pos, one too small
        ("31ZP", (0,0,0,0)),  # numbers too small
        ("02NE", (0,0,0,0)),  # numbers too small
        ("36ZF", (1,0,0,0)),  # one number correct pos, one too small
        ("26SR", (1,0,0,1)),  # one number correct pos, one small, one letter wrong pos
        ("65FG", (0,1,0,1))   # one number wrong pos, one small, one letter correct pos
    ]
    
    for guess, expected in guesses_and_feedback:
        result = verify_guess(password, guess, None)
        if result != expected:
            return False
    return True

# We know the numbers are 7,6 from condition 6
password = ['7', '6', 'R', 'G']
print(f"Testing password {password}")
print(f"Valid: {test_password(password)}")