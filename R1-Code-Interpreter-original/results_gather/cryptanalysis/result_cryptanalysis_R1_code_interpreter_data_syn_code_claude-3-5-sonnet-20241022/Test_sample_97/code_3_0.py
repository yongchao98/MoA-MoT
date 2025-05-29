def analyze_conditions():
    # Initialize possible values
    possible_numbers = set(range(0, 10))
    possible_letters = set([chr(i) for i in range(ord('A'), ord('Z')+1)])
    
    # Analyze each guess and its feedback
    
    # From "71BQ" - one number correct but wrong pos, one small, both letters early
    # This means one of 7,1 is in password but wrong pos, other number > other digit
    # B,Q are both before correct letters
    possible_first_nums = set()
    possible_second_nums = set()
    
    # From "96MJ" - one number correct in correct pos, one small
    # This means 9 or 6 is in correct position
    # From multiple conditions about numbers being too small
    # We know larger numbers are needed
    possible_first_nums.add(9)
    possible_second_nums.add(8)
    
    # From "40QT" - both numbers too small, Q or T correct but wrong pos
    # This confirms Q,T are in the password
    possible_letters = possible_letters.intersection({'Q', 'T'})
    
    # From multiple "both letters early" feedback
    # Letters must be late in alphabet
    
    # From "89XJ" - one number correct but wrong pos
    # This confirms 8 is in the password
    
    # Combine all conditions
    passwords = []
    for n1 in [9]:  # 9 must be first based on "96MJ"
        for n2 in [8]:  # 8 must be second based on multiple conditions
            for l1 in ['Q', 'T']:
                for l2 in ['Q', 'T']:
                    if l1 != l2:
                        password = [str(n1), str(n2), l1, l2]
                        # Verify against key conditions
                        if verify_conditions(password):
                            passwords.append(password)
    return passwords

def verify_conditions(password):
    # Helper function to check if a guess matches conditions
    def check_guess(guess, password, num_correct_pos, num_wrong_pos, num_small, letters_early):
        g_nums = [int(guess[0]), int(guess[1])]
        g_letters = [guess[2], guess[3]]
        p_nums = [int(password[0]), int(password[1])]
        p_letters = [password[2], password[3]]
        
        # Check numbers in correct position
        correct_pos = sum(1 for i in range(2) if g_nums[i] == p_nums[i])
        if correct_pos != num_correct_pos:
            return False
            
        # Check numbers in wrong position
        wrong_pos = sum(1 for i in range(2) for j in range(2) 
                       if i != j and g_nums[i] == p_nums[j])
        if wrong_pos != num_wrong_pos:
            return False
            
        # Check small numbers
        small = sum(1 for i in range(2) 
                   if g_nums[i] < p_nums[i] and g_nums[i] not in p_nums)
        if small != num_small:
            return False
            
        # Check early letters
        if letters_early:
            if not all(g_letters[i] < p_letters[i] for i in range(2)):
                return False
                
        return True
    
    # Check key guesses
    if not check_guess("96MJ", password, 1, 0, 1, True):  # One correct pos, one small
        return False
    if not check_guess("40QT", password, 0, 0, 2, False):  # Both small, Q/T wrong pos
        return False
    if not check_guess("89XJ", password, 0, 1, 0, False):  # One correct wrong pos
        return False
    
    return True

result = analyze_conditions()
print(result)