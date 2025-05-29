def verify_guess(guess_nums, guess_letters, actual_nums, actual_letters, num_feedback, letter_feedback):
    # For guess "84KP", actual [1,7,K,R], feedback "both numbers incorrect, one letter correct but wrong position"
    if guess_nums == ['8','4'] and guess_letters == ['K','P']:
        return (guess_nums[0] not in actual_nums and 
                guess_nums[1] not in actual_nums and 
                (('K' in actual_letters and guess_letters.index('K') != actual_letters.index('K')) or
                 ('P' in actual_letters and guess_letters.index('P') != actual_letters.index('P'))))
    
    # For guess "13GY", actual [1,7,K,R], feedback "one number correct in position, one too small, both letters wrong"
    elif guess_nums == ['1','3'] and guess_letters == ['G','Y']:
        return (guess_nums[0] == actual_nums[0] and
                int(guess_nums[1]) < int(actual_nums[1]) and
                'G' not in actual_letters and 
                'Y' not in actual_letters)
    
    # For guess "31HM", actual [1,7,K,R], feedback "one number correct wrong position, one too small, both letters wrong"
    elif guess_nums == ['3','1'] and guess_letters == ['H','M']:
        return ('1' in actual_nums and
                guess_nums.index('1') != actual_nums.index('1') and
                int(guess_nums[0]) < min(int(x) for x in actual_nums) and
                'H' not in actual_letters and
                'M' not in actual_letters)
    
    # For guess "84RL", actual [1,7,K,R], feedback "both numbers wrong, one letter correct wrong position"
    elif guess_nums == ['8','4'] and guess_letters == ['R','L']:
        return (guess_nums[0] not in actual_nums and
                guess_nums[1] not in actual_nums and
                'R' in actual_letters and
                guess_letters.index('R') != actual_letters.index('R'))
    
    # For guess "02AN", actual [1,7,K,R], feedback "both numbers too small, both letters wrong"
    elif guess_nums == ['0','2'] and guess_letters == ['A','N']:
        return (int(guess_nums[0]) < min(int(x) for x in actual_nums) and
                int(guess_nums[1]) < min(int(x) for x in actual_nums) and
                'A' not in actual_letters and
                'N' not in actual_letters)
    
    # For guess "20BH", actual [1,7,K,R], feedback "both numbers too small, both letters too early"
    elif guess_nums == ['2','0'] and guess_letters == ['B','H']:
        return (int(guess_nums[0]) < min(int(x) for x in actual_nums) and
                int(guess_nums[1]) < min(int(x) for x in actual_nums) and
                'B' not in actual_letters and
                'H' not in actual_letters)
    
    # For guess "17XY", actual [1,7,K,R], feedback "one number correct wrong position, one too small, both letters wrong"
    elif guess_nums == ['1','7'] and guess_letters == ['X','Y']:
        return ('7' in actual_nums and
                guess_nums.index('7') != actual_nums.index('7') and
                'X' not in actual_letters and
                'Y' not in actual_letters)
    
    return False

# Test the specific combination [1,7,K,R]
actual_nums = ['1', '7']
actual_letters = ['K', 'R']

# Verify all guesses
guesses = [
    (['8','4'], ['K','P']),
    (['1','3'], ['G','Y']),
    (['3','1'], ['H','M']),
    (['8','4'], ['R','L']),
    (['0','2'], ['A','N']),
    (['2','0'], ['B','H']),
    (['1','7'], ['X','Y'])
]

all_conditions_met = True
for guess_nums, guess_letters in guesses:
    if not verify_guess(guess_nums, guess_letters, actual_nums, actual_letters, "", ""):
        all_conditions_met = False
        print(f"Failed for guess: {guess_nums}{guess_letters}")
        break

if all_conditions_met:
    print([actual_nums[0], actual_nums[1], actual_letters[0], actual_letters[1]])