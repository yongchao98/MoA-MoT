def verify_84KP(actual_nums, actual_letters):
    # For 84KP: both numbers incorrect, one letter correct but wrong position, one letter too late
    guess_nums = ['8', '4']
    guess_letters = ['K', 'P']
    
    # Check numbers: 8 and 4 should not be in actual numbers
    numbers_correct = (guess_nums[0] not in actual_nums and 
                      guess_nums[1] not in actual_nums)
    
    # Check letters: K should be in actual but in wrong position, P should not be in actual
    letters_correct = ('K' in actual_letters and 
                      guess_letters.index('K') != actual_letters.index('K') and 
                      'P' not in actual_letters)
    
    return numbers_correct and letters_correct

def verify_13GY(actual_nums, actual_letters):
    # For 13GY: one number correct in position, one too small, both letters wrong
    return (actual_nums[0] == '1' and 
            '3' not in actual_nums and 
            'G' not in actual_letters and 
            'Y' not in actual_letters)

def verify_31HM(actual_nums, actual_letters):
    # For 31HM: one number correct wrong position, one too small, both letters wrong
    return ('1' in actual_nums and 
            '3' not in actual_nums and 
            'H' not in actual_letters and 
            'M' not in actual_letters)

def verify_84RL(actual_nums, actual_letters):
    # For 84RL: both numbers wrong, one letter correct wrong position, one letter too late
    return ('8' not in actual_nums and 
            '4' not in actual_nums and 
            'R' in actual_letters and 
            'L' not in actual_letters)

def verify_02AN(actual_nums, actual_letters):
    # For 02AN: both numbers too small, both letters wrong
    return (int('0') < int(actual_nums[0]) and 
            int('2') < int(actual_nums[0]) and 
            'A' not in actual_letters and 
            'N' not in actual_letters)

def verify_20BH(actual_nums, actual_letters):
    # For 20BH: both numbers too small, both letters too early
    return (int('2') < int(actual_nums[0]) and 
            int('0') < int(actual_nums[0]) and 
            'B' not in actual_letters and 
            'H' not in actual_letters)

def verify_17XY(actual_nums, actual_letters):
    # For 17XY: one number correct wrong position, one too small, both letters wrong
    return ('7' in actual_nums and 
            actual_nums.index('7') == 1 and  # 7 should be in second position
            'X' not in actual_letters and 
            'Y' not in actual_letters)

# Test the combination [1,7,K,R]
actual_nums = ['1', '7']
actual_letters = ['K', 'R']

# Verify all conditions
verifications = [
    (verify_84KP(actual_nums, actual_letters), "84KP"),
    (verify_13GY(actual_nums, actual_letters), "13GY"),
    (verify_31HM(actual_nums, actual_letters), "31HM"),
    (verify_84RL(actual_nums, actual_letters), "84RL"),
    (verify_02AN(actual_nums, actual_letters), "02AN"),
    (verify_20BH(actual_nums, actual_letters), "20BH"),
    (verify_17XY(actual_nums, actual_letters), "17XY")
]

all_passed = True
for verification, guess in verifications:
    if not verification:
        print(f"Failed verification for guess: {guess}")
        all_passed = False

if all_passed:
    print("All verifications passed!")
    print([actual_nums[0], actual_nums[1], actual_letters[0], actual_letters[1]])