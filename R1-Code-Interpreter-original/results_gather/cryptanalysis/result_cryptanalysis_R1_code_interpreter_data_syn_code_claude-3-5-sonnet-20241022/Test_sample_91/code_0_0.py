def check_guess(guess_nums, guess_letters, conditions):
    # Check if numbers or letters repeat
    if len(set(guess_nums)) != 2 or len(set(guess_letters)) != 2:
        return False
    
    # Check condition 1: 47LE
    if guess_nums[0] == 4 or guess_nums[1] == 7:
        return False
    
    # Check condition 2: 10FX
    correct_wrong_pos = sum(n in [1,0] for n in guess_nums)
    if correct_wrong_pos != 1:
        return False
    
    # Check condition 3: 50ZK
    correct_wrong_pos = sum(n in [5,0] for n in guess_nums)
    if correct_wrong_pos != 1:
        return False
    
    # Check condition 4: 48IM
    if not (set(guess_letters) == set(['I', 'M'])):
        return False
    
    # Check condition 5: 86NM
    if 'M' not in guess_letters:
        return False
    
    # Check condition 6: 39KL
    correct_wrong_pos = sum(n in [3,9] for n in guess_nums)
    if correct_wrong_pos != 1:
        return False
    
    # Check all numbers are less than 4 (from multiple conditions)
    if any(n >= 4 for n in guess_nums):
        return False
        
    return True

# Try all possible combinations
for n1 in range(4):
    for n2 in range(4):
        if n1 == n2:
            continue
        for l1 in ['I', 'M']:
            for l2 in ['I', 'M']:
                if l1 == l2:
                    continue
                guess = [n1, n2, l1, l2]
                if check_guess(guess[:2], guess[2:], None):
                    print(guess)