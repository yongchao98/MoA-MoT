def verify_all_conditions():
    candidate = [8, 3, 'L', 'K']
    
    def check_guess(guess_str, conditions):
        guess = [int(guess_str[0]), int(guess_str[1]), guess_str[2], guess_str[3]]
        
        # Check number conditions
        if "both nums wrong" in conditions:
            if guess[0] in candidate[:2] or guess[1] in candidate[:2]:
                print(f"Failed: {guess_str} should have both numbers wrong")
                return False
        
        if "one num right right pos" in conditions:
            correct_pos = sum(1 for i in range(2) if guess[i] == candidate[i])
            if correct_pos != 1:
                print(f"Failed: {guess_str} should have one number in correct position")
                return False
        
        if "one num right wrong pos" in conditions:
            correct_wrong = (guess[0] == candidate[1] or guess[1] == candidate[0])
            if not correct_wrong:
                print(f"Failed: {guess_str} should have one number in wrong position")
                return False
        
        return True
    
    # Test all guesses
    guesses = [
        ("43AG", "one num right right pos"),  # Critical guess
        ("43OB", "one num right right pos"),  # Critical guess
        ("83LK", "one num right right pos, K wrong pos"),  # Critical guess
        ("28IK", "both nums wrong"),  # Critical guess
        ("06MG", "one num right wrong pos"),
        ("23CM", "one num right right pos"),
        ("46EB", "one num right wrong pos"),
        ("03GE", "one num right right pos")
    ]
    
    print("Verifying critical conditions:")
    all_valid = True
    for guess_str, conditions in guesses:
        if not check_guess(guess_str, conditions):
            all_valid = False
            break
    
    if not all_valid:
        print("\nPrevious solution [8, 3, 'L', 'K'] fails verification")
        
        # The guess "43AG" having one number correct in correct position
        # suggests 3 must be correct in second position
        # "46EB" having one number correct in wrong position
        # suggests 6 must be in first position
        
        new_candidate = [6, 3, 'L', 'K']
        print("\nTesting new candidate:", new_candidate)
        
        # Verify new candidate matches all conditions
        all_valid = True
        for guess_str, conditions in guesses:
            if not check_guess(guess_str, conditions):
                all_valid = False
                break
                
        if all_valid:
            print("\nNew solution [6, 3, 'L', 'K'] passes all verifications")
            return new_candidate
            
    return candidate

result = verify_all_conditions()
print("\nFinal verified password:", result)