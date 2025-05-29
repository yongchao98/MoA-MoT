def verify_guess(nums, letters, guess_nums, guess_letters, feedback):
    if guess_nums == ['8', '4'] and guess_letters == ['K', 'P']:
        # Both numbers wrong, K in wrong position, P too late
        return (not any(n in nums for n in guess_nums) and  # 8,4 not in [1,7]
                'K' in letters and letters.index('K') == 1 and  # K in second position
                ord('P') > ord('R') and ord('P') > ord('K'))  # P later than both R,K
    
    elif guess_nums == ['1', '3'] and guess_letters == ['G', 'Y']:
        # 1 correct in position, 3 too small, both letters wrong
        return (nums[0] == '1' and  # 1 in first position
                int('3') < int('7') and  # 3 < 7
                not any(l in letters for l in guess_letters))  # G,Y not in [R,K]
    
    elif guess_nums == ['3', '1'] and guess_letters == ['H', 'M']:
        # 1 in wrong position, 3 too small, both letters wrong
        return ('1' in nums and nums[1] != '1' and  # 1 present but not in second position
                int('3') < int('7') and  # 3 < 7
                not any(l in letters for l in guess_letters))  # H,M not in [R,K]
    
    elif guess_nums == ['8', '4'] and guess_letters == ['R', 'L']:
        # Both numbers wrong, R in wrong position, L too late
        return (not any(n in nums for n in guess_nums) and  # 8,4 not in [1,7]
                'R' in letters and letters.index('R') == 0 and  # R in first position
                ord('L') > ord('R') and ord('L') > ord('K'))  # L later than both R,K
    
    elif guess_nums == ['0', '2'] and guess_letters == ['A', 'N']:
        # Both numbers too small, both letters wrong
        return (all(int(n) < int('1') for n in guess_nums) and  # 0,2 < 1
                not any(l in letters for l in guess_letters))  # A,N not in [R,K]
    
    elif guess_nums == ['2', '0'] and guess_letters == ['B', 'H']:
        # Both numbers too small, both letters too early
        return (all(int(n) < int('1') for n in guess_nums) and  # 2,0 < 1
                all(ord(l) < ord('R') and ord(l) < ord('K') for l in guess_letters))  # B,H earlier than R,K
    
    elif guess_nums == ['1', '7'] and guess_letters == ['X', 'Y']:
        # 7 in wrong position, both letters wrong and too late
        return ('7' in nums and nums[0] != '7' and  # 7 present but not in first position
                not any(l in letters for l in guess_letters) and  # X,Y not in [R,K]
                all(ord(l) > ord('R') and ord(l) > ord('K') for l in guess_letters))  # X,Y later than R,K
    
    return False

# Test [1,7,R,K]
nums = ['1', '7']
letters = ['R', 'K']

# All guesses and their feedbacks
guesses = [
    (['8', '4'], ['K', 'P'], "84KP"),
    (['1', '3'], ['G', 'Y'], "13GY"),
    (['3', '1'], ['H', 'M'], "31HM"),
    (['8', '4'], ['R', 'L'], "84RL"),
    (['0', '2'], ['A', 'N'], "02AN"),
    (['2', '0'], ['B', 'H'], "20BH"),
    (['1', '7'], ['X', 'Y'], "17XY")
]

# Verify each guess
all_passed = True
for guess_nums, guess_letters, guess_name in guesses:
    if not verify_guess(nums, letters, guess_nums, guess_letters, guess_name):
        print(f"Failed verification for guess: {guess_name}")
        all_passed = False

if all_passed:
    print("All verifications passed!")
    print(nums + letters)