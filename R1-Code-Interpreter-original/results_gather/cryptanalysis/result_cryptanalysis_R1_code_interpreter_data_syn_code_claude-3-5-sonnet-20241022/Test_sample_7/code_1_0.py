def check_guess(candidate, guess, conditions):
    # Split into numbers and letters
    cand_nums = candidate[:2]
    cand_lets = candidate[2:]
    guess_nums = guess[:2]
    guess_lets = guess[2:]
    
    # Count conditions
    num_correct_pos = 0
    num_too_large = 0
    let_correct_pos = 0
    let_too_late = 0
    let_wrong = 0
    
    # Check numbers
    for i in range(2):
        if guess_nums[i] == cand_nums[i]:
            num_correct_pos += 1
        elif int(guess_nums[i]) > int(cand_nums[i]):
            num_too_large += 1
            
    # Check letters
    for i in range(2):
        if guess_lets[i] == cand_lets[i]:
            let_correct_pos += 1
        elif guess_lets[i] > cand_lets[i]:
            let_too_late += 1
        else:
            let_wrong += 1
            
    return (num_correct_pos, num_too_large, let_correct_pos, let_too_late, let_wrong)

def find_password():
    # Define the guesses and their conditions
    guesses = [
        ('08VY', (1, 1, 0, 2, 0)),  # one num correct pos, one too large, both letters too late
        ('07CX', (1, 1, 0, 0, 2)),  # one num correct pos, one too large, both letters wrong
        ('45NJ', (0, 2, 1, 0, 1)),  # both nums too large, one letter correct pos, one letter wrong
        ('28YQ', (0, 1, 0, 1, 1))   # one num wrong pos(not counted here), one too large, one letter wrong pos(not counted), one too late
    ]
    
    import itertools
    
    # Generate all possible combinations
    numbers = '0123456789'
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    for nums in itertools.permutations(numbers, 2):
        for lets in itertools.permutations(letters, 2):
            candidate = list(nums + lets)
            
            # Check if this candidate matches all guesses
            matches_all = True
            for guess, conditions in guesses:
                if check_guess(candidate, guess, conditions) != conditions:
                    matches_all = False
                    break
                    
                # Additional checks for guess 4 (28YQ)
                if guess == '28YQ':
                    # Check if 2 is in candidate numbers but wrong position
                    if '2' not in nums or '2' == candidate[nums.index('2')]:
                        matches_all = False
                        break
                    # Check if Y is in candidate letters but wrong position
                    if 'Y' not in lets or 'Y' == candidate[2 + lets.index('Y')]:
                        matches_all = False
                        break
            
            if matches_all:
                print(candidate)
                return candidate

result = find_password()