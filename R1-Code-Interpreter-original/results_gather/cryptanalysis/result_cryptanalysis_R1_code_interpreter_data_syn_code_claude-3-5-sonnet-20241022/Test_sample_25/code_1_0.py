def verify_solution(solution, verbose=True):
    def check_guess(guess, feedback, solution):
        # Split into numbers and letters
        guess_nums = [int(guess[0]), int(guess[1])]
        guess_letters = [guess[2], guess[3]]
        sol_nums = [int(solution[0]), int(solution[1])]
        sol_letters = [solution[2], solution[3]]
        
        # Count correct numbers
        nums_correct_pos = sum(1 for i in range(2) if guess_nums[i] == sol_nums[i])
        nums_correct_val = sum(1 for x in guess_nums if x in sol_nums) - nums_correct_pos
        nums_too_small = sum(1 for i in range(2) if guess_nums[i] < sol_nums[i])
        
        # Count correct letters
        letters_correct_pos = sum(1 for i in range(2) if guess_letters[i] == sol_letters[i])
        letters_correct_val = sum(1 for x in guess_letters if x in sol_letters) - letters_correct_pos
        
        if verbose:
            print(f"\nGuess: {guess}")
            print(f"Numbers: correct pos = {nums_correct_pos}, correct val = {nums_correct_val}, too small = {nums_too_small}")
            print(f"Letters: correct pos = {letters_correct_pos}, correct val = {letters_correct_val}")
        
        # Check against feedback
        if feedback == 1:  # 91SF
            return (nums_correct_pos == 0 and 
                   letters_correct_pos == 1 and letters_correct_val == 0)
        elif feedback == 2:  # 03HF
            return (nums_correct_pos == 1 and nums_too_small == 1 and 
                   letters_correct_pos == 0 and letters_correct_val == 0)
        elif feedback == 3:  # 75CP
            return (nums_correct_pos == 0 and 
                   letters_correct_pos == 0 and letters_correct_val == 1)
        elif feedback == 4:  # 64YB
            return (nums_correct_pos == 1 and nums_correct_val == 0 and 
                   letters_correct_pos == 0 and letters_correct_val == 0)
        
        return False

    guesses = [
        ('91SF', 1),
        ('03HF', 2),
        ('75CP', 3),
        ('64YB', 4)
    ]
    
    print(f"\nVerifying solution: {solution}")
    all_match = True
    for guess, feedback in guesses:
        matches = check_guess(guess, feedback, solution)
        if verbose:
            print(f"Matches feedback {feedback}: {matches}")
        if not matches:
            all_match = False
    
    return all_match

# Test both solutions
solutions = [['0', '4', 'S', 'C'], ['6', '3', 'S', 'C']]
for solution in solutions:
    verify_solution(solution)