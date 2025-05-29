def is_valid_password(password, guesses_feedback):
    # Helper function to check if a letter is "too early"
    def is_too_early(guess_letter, actual_letter):
        return guess_letter < actual_letter
    
    # Check each guess against the password
    for guess, feedback in guesses_feedback:
        # Split numbers and letters
        guess_nums = [int(guess[0]), int(guess[1])]
        guess_letters = [guess[2], guess[3]]
        pass_nums = [int(password[0]), int(password[1])]
        pass_letters = [password[2], password[3]]
        
        if feedback == "one_correct_wrong_pos_num_one_small_both_letters_early":
            # Check numbers
            correct_wrong_pos = False
            one_small = False
            for i, gn in enumerate(guess_nums):
                for j, pn in enumerate(pass_nums):
                    if gn == pn and i != j:
                        correct_wrong_pos = True
                    elif gn < pn:
                        one_small = True
            if not (correct_wrong_pos and one_small):
                return False
            # Check letters
            if not (guess_letters[0] < pass_letters[0] and guess_letters[1] < pass_letters[1]):
                return False
                
        elif feedback == "both_nums_small_one_letter_correct_wrong_pos_one_early":
            # Check numbers
            if not (guess_nums[0] < pass_nums[0] and guess_nums[1] < pass_nums[1]):
                return False
            # Check letters
            correct_wrong_pos = False
            for i, gl in enumerate(guess_letters):
                for j, pl in enumerate(pass_letters):
                    if gl == pl and i != j:
                        correct_wrong_pos = True
            if not correct_wrong_pos:
                return False
                
        elif feedback == "one_num_correct_pos_one_small_both_letters_early":
            # Check numbers
            correct_pos = False
            one_small = False
            for i in range(2):
                if guess_nums[i] == pass_nums[i]:
                    correct_pos = True
                elif guess_nums[i] < pass_nums[i]:
                    one_small = True
            if not (correct_pos and one_small):
                return False
            # Check letters
            if not (guess_letters[0] < pass_letters[0] and guess_letters[1] < pass_letters[1]):
                return False
    
    return True

# Test all possible combinations
numbers = range(10)
letters = [chr(i) for i in range(65, 91)]  # A-Z
valid_passwords = []

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # Numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # Letters can't repeat
                        password = [str(n1), str(n2), l1, l2]
                        # Check against known conditions
                        if (int(password[0]) > 7 or int(password[1]) > 7) and \
                           l1 > 'M' and l2 > 'M':  # Based on various conditions
                            valid_passwords.append(password)

print(valid_passwords)