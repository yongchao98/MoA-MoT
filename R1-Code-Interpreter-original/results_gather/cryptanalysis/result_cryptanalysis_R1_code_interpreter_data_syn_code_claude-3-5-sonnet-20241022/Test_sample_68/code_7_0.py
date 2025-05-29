def verify_84KP(actual_nums, actual_letters):
    # For 84KP: 
    # - both numbers incorrect
    numbers_ok = '8' not in actual_nums and '4' not in actual_nums
    
    # - one letter (K) correct but wrong position
    k_ok = 'K' in actual_letters and actual_letters.index('K') != 0  # K should be present but not in same position
    
    # - P is too late in alphabet
    p_ok = ord('P') > ord(actual_letters[0]) and ord('P') > ord(actual_letters[1])
    
    return numbers_ok and k_ok and p_ok

def verify_13GY(actual_nums, actual_letters):
    # - 1 correct in position
    # - 3 too small
    # - both letters wrong
    return (actual_nums[0] == '1' and 
            int('3') < int(actual_nums[1]) and 
            'G' not in actual_letters and 
            'Y' not in actual_letters)

def verify_31HM(actual_nums, actual_letters):
    # - 1 correct but wrong position
    # - 3 too small
    # - both letters wrong
    return ('1' in actual_nums and 
            actual_nums.index('1') != 1 and
            int('3') < int(actual_nums[1]) and 
            'H' not in actual_letters and 
            'M' not in actual_letters)

def verify_84RL(actual_nums, actual_letters):
    # - both numbers incorrect
    numbers_ok = '8' not in actual_nums and '4' not in actual_nums
    
    # - one letter (R) correct but wrong position
    r_ok = 'R' in actual_letters and actual_letters.index('R') != 0
    
    # - L is too late in alphabet
    l_ok = ord('L') > ord(actual_letters[0]) and ord('L') > ord(actual_letters[1])
    
    return numbers_ok and r_ok and l_ok

def verify_02AN(actual_nums, actual_letters):
    # - both numbers too small
    numbers_ok = all(int(n) < min(int(actual_nums[0]), int(actual_nums[1])) for n in ['0', '2'])
    
    # - both letters wrong
    letters_ok = 'A' not in actual_letters and 'N' not in actual_letters
    
    return numbers_ok and letters_ok

def verify_20BH(actual_nums, actual_letters):
    # - both numbers too small
    numbers_ok = all(int(n) < min(int(actual_nums[0]), int(actual_nums[1])) for n in ['2', '0'])
    
    # - both letters too early
    letters_ok = (ord('B') < min(ord(actual_letters[0]), ord(actual_letters[1])) and 
                 ord('H') < min(ord(actual_letters[0]), ord(actual_letters[1])))
    
    return numbers_ok and letters_ok

def verify_17XY(actual_nums, actual_letters):
    # - 7 correct but wrong position
    seven_ok = '7' in actual_nums and actual_nums.index('7') != 0
    
    # - 1 wrong and too small
    one_ok = int('1') < int(actual_nums[1])
    
    # - both letters wrong and too late
    letters_ok = ('X' not in actual_letters and 'Y' not in actual_letters and
                 all(ord(l) > max(ord(actual_letters[0]), ord(actual_letters[1])) for l in ['X', 'Y']))
    
    return seven_ok and one_ok and letters_ok

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

# Test the combination
all_passed = True
for verification, guess in verifications:
    if not verification:
        print(f"Failed verification for guess: {guess}")
        all_passed = False

if all_passed:
    print("All verifications passed!")
    print([actual_nums[0], actual_nums[1], actual_letters[0], actual_letters[1]])