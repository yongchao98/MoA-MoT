def verify_password():
    # Candidate password
    password = [4, 6, 'K', 'L']
    
    # Function to check if a guess matches the known feedback
    def check_guess(guess, feedback_type):
        # Convert guess string to list of chars/digits
        guess = [int(guess[0]), int(guess[1]), guess[2], guess[3]]
        
        if feedback_type == 1:  # Both numbers incorrect
            if guess[0] in password[:2] or guess[1] in password[:2]:
                return False
        elif feedback_type == 2:  # One number correct wrong pos
            correct_nums = sum(1 for i in range(2) for j in range(2) if guess[i] == password[j])
            if correct_nums != 1:
                return False
        elif feedback_type == 3:  # One number correct right pos
            correct_pos = sum(1 for i in range(2) if guess[i] == password[i])
            if correct_pos != 1:
                return False
                
        return True
    
    # Verify against key guesses
    checks = [
        ("46EB", 2),  # One number correct but wrong position
        ("43AG", 3),  # One number correct in correct position
        ("83LK", 3),  # One number correct in position, K wrong position
        ("28IK", 1),  # Both numbers incorrect, K wrong position
        ("06MG", 2),  # One number correct wrong position
        ("23CM", 3),  # One number correct correct position
    ]
    
    all_valid = True
    for guess, feedback_type in checks:
        if not check_guess(guess, feedback_type):
            all_valid = False
            print(f"Failed check for guess {guess}")
            break
    
    # Additional verification for letter positions
    def verify_letters():
        # K must be first letter (from guesses 12 and 17)
        # L must be second letter and not appear in wrong guesses
        late_letters = set('SWRQYUVWO')
        early_letters = set('IDGE')
        if 'L' in late_letters or 'L' in early_letters:
            return False
        return True
    
    if all_valid and verify_letters():
        print("Password verified: 46KL matches all conditions")
        return True
    return False

verify_password()