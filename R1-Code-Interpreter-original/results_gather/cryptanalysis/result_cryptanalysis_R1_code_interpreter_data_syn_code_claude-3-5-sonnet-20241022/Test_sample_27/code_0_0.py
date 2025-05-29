def verify_password(password, guesses):
    # Helper function to check if a guess matches the feedback
    def check_guess(guess, feedback):
        if feedback == "both_numbers_incorrect":
            return (password[0] != int(guess[0]) and 
                   password[1] != int(guess[1]))
        elif feedback == "both_numbers_correct":
            return (password[0] == int(guess[0]) and 
                   password[1] == int(guess[1]))
        elif feedback == "numbers_too_small":
            return (password[0] > int(guess[0]) and 
                   password[1] > int(guess[1]))
        elif feedback == "one_number_correct_one_small":
            return ((password[0] == int(guess[0]) and int(guess[1]) < password[1]) or
                   (password[1] == int(guess[1]) and int(guess[0]) < password[0]))
        elif feedback == "letters_incorrect":
            return (password[2] != guess[2] and 
                   password[3] != guess[3])
        elif feedback == "letters_too_early":
            return (password[2] > guess[2] and 
                   password[3] > guess[3])
        elif feedback == "letters_correct_wrong_pos":
            letters = {guess[2], guess[3]}
            return (set([password[2], password[3]]) == letters and
                   (password[2] != guess[2] or password[3] != guess[3]))

    # Check all conditions
    conditions = [
        ("92OR", "both_numbers_incorrect", "letters_incorrect"),
        ("86YQ", "both_numbers_correct", "letters_incorrect"),
        ("20CF", "numbers_too_small", "letters_too_early"),
        ("80ZI", "one_number_correct_one_small", "letters_correct_wrong_pos")
    ]
    
    for guess, num_feedback, letter_feedback in conditions:
        if not (check_guess(guess, num_feedback) and 
                check_guess(guess, letter_feedback)):
            return False
    return True

# Test all possible combinations
numbers = [8, 6]  # We know these from guess 2
letters = ['Z', 'I']  # We know these from guess 4

# Try both possible arrangements of letters
possible_passwords = []
for l1, l2 in [(letters[0], letters[1]), (letters[1], letters[0])]:
    password = [8, 6, l1, l2]
    if verify_password(password, []):
        possible_passwords.append(password)

print(possible_passwords)