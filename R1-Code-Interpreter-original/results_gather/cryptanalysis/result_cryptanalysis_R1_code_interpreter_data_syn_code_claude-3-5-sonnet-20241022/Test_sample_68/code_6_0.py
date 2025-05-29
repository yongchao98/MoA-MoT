def verify_84KP(actual_nums, actual_letters):
    # For 84KP: both numbers incorrect, one letter correct but wrong position, one letter too late
    # Numbers: 8,4 not in [1,7]
    numbers_correct = '8' not in actual_nums and '4' not in actual_nums
    
    # Letters: K in [K,R] but wrong position, P too late in alphabet than both K and R
    letters_correct = ('K' in actual_letters and 
                      ord('P') > ord('K') and ord('P') > ord('R'))
    
    return numbers_correct and letters_correct

def verify_13GY(actual_nums, actual_letters):
    # For 13GY: one number correct in position, one too small, both letters wrong
    return (actual_nums[0] == '1' and 
            int('3') < int(actual_nums[1]) and 
            'G' not in actual_letters and 
            'Y' not in actual_letters)

def verify_31HM(actual_nums, actual_letters):
    # For 31HM: one number correct wrong position, one too small, both letters wrong
    return ('1' in actual_nums and 
            int('3') < int(actual_nums[1]) and 
            'H' not in actual_letters and 
            'M' not in actual_letters)

def verify_84RL(actual_nums, actual_letters):
    # For 84RL: both numbers wrong, one letter correct wrong position, one letter too late
    return ('8' not in actual_nums and 
            '4' not in actual_nums and 
            'R' in actual_letters and 
            ord('L') > ord('K') and ord('L') > ord('R'))

def verify_02AN(actual_nums, actual_letters):
    # For 02AN: both numbers too small, both letters wrong
    min_actual = min(int(actual_nums[0]), int(actual_nums[1]))
    return (int('0') < min_actual and 
            int('2') < min_actual and 
            'A' not in actual_letters and 
            'N' not in actual_letters)

def verify_20BH(actual_nums, actual_letters):
    # For 20BH: both numbers too small, both letters too early
    min_actual = min(int(actual_nums[0]), int(actual_nums[1]))
    min_letter = min(ord('K'), ord('R'))
    return (int('2') < min_actual and 
            int('0') < min_actual and 
            ord('B') < min_letter and 
            ord('H') < min_letter)

def verify_17XY(actual_nums, actual_letters):
    # For 17XY: one number correct wrong position, one too small, both letters wrong and too late
    max_letter = max(ord('K'), ord('R'))
    return ('7' in actual_nums and 
            actual_nums.index('7') == 1 and
            'X' not in actual_letters and 
            'Y' not in actual_letters and
            ord('X') > max_letter and 
            ord('Y') > max_letter)

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