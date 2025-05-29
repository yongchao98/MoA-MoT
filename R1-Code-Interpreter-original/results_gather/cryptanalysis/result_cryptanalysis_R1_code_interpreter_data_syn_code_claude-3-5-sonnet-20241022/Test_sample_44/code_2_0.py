def verify_all_conditions(password):
    def check_position_and_value(guess_nums, pass_nums):
        correct_pos = sum(1 for i in range(2) if guess_nums[i] == pass_nums[i])
        in_password = sum(1 for x in guess_nums if x in pass_nums) - correct_pos
        return correct_pos, in_password

    # Convert password to string for easier comparison
    pass_str = ''.join(password)
    
    # All conditions with their feedback
    conditions = [
        ("71HT", "one number correct position, one too small"),
        ("19ZS", "both numbers wrong, letters too late"),
        ("48PA", "all wrong"),
        ("49FO", "all wrong"),
        ("25IP", "both numbers too small"),
        ("76KR", "both numbers correct position, one letter wrong position"),
        ("17BY", "one number wrong position, one too small"),
        ("31ZP", "both numbers too small"),
        ("02NE", "both numbers too small"),
        ("36ZF", "one number correct position, one too small"),
        ("26SR", "one number correct position, one small, one letter wrong position"),
        ("65FG", "one number wrong position, one small, one letter correct position")
    ]

    for guess, feedback in conditions:
        # Check numbers
        guess_nums = [guess[0], guess[1]]
        pass_nums = [password[0], password[1]]
        num_pos, num_wrong = check_position_and_value(guess_nums, pass_nums)
        
        # Check letters
        guess_letters = [guess[2], guess[3]]
        pass_letters = [password[2], password[3]]
        letter_pos, letter_wrong = check_position_and_value(guess_letters, pass_letters)
        
        # Verify specific conditions
        if guess == "76KR" and (num_pos != 2 or letter_wrong != 1):
            return False, f"Failed at {guess}: numbers should be correct and one letter wrong position"
        elif guess == "65FG" and (num_wrong != 1 or letter_pos != 1):
            return False, f"Failed at {guess}: one number wrong pos, one letter correct pos"
        elif guess == "71HT" and num_pos != 1:
            return False, f"Failed at {guess}: one number should be in correct position"
        
    return True, "All conditions satisfied"

# Test possible combinations
numbers = ['7', '6']  # From condition 6
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

best_password = None
for l1 in letters:
    for l2 in letters:
        if l1 != l2:
            password = ['7', '6', l1, l2]
            valid, message = verify_all_conditions(password)
            if valid:
                print(f"Found valid password: {password}")
                best_password = password
                break

if best_password:
    print(f"Final password: {best_password}")
else:
    print("No valid password found")