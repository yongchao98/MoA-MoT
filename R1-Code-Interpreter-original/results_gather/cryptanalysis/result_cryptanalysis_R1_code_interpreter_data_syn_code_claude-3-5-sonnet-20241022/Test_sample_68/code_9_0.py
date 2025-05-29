def verify_all_conditions(nums, letters):
    # Helper functions
    def is_too_late(letter, password_letters):
        return ord(letter) > ord(password_letters[0]) and ord(letter) > ord(password_letters[1])
    
    def is_too_early(letter, password_letters):
        return ord(letter) < ord(password_letters[0]) and ord(letter) < ord(password_letters[1])
    
    # Verify each guess individually with exact conditions
    
    # Guess 1: 84KP
    if not (
        '8' not in nums and '4' not in nums and  # both numbers wrong
        'K' in letters and  # K must be in password
        'P' not in letters and  # P not in password
        is_too_late('P', letters)  # P must be too late
    ):
        return "Failed 84KP"
    
    # Guess 2: 13GY
    if not (
        nums[0] == '1' and  # 1 correct in first position
        '3' not in nums and int('3') < int(nums[1]) and  # 3 too small
        'G' not in letters and 'Y' not in letters  # both letters wrong
    ):
        return "Failed 13GY"
    
    # Guess 3: 31HM
    if not (
        '1' in nums and nums[1] != '1' and  # 1 in wrong position
        '3' not in nums and int('3') < int(nums[1]) and  # 3 too small
        'H' not in letters and 'M' not in letters  # both letters wrong
    ):
        return "Failed 31HM"
    
    # Guess 4: 84RL
    if not (
        '8' not in nums and '4' not in nums and  # both numbers wrong
        'R' in letters and  # R must be in password
        'L' not in letters and  # L not in password
        is_too_late('L', letters)  # L must be too late
    ):
        return "Failed 84RL"
    
    # Guess 5: 02AN
    if not (
        all(int(n) < int(nums[0]) and int(n) < int(nums[1]) for n in ['0', '2']) and  # both numbers too small
        'A' not in letters and 'N' not in letters  # both letters wrong
    ):
        return "Failed 02AN"
    
    # Guess 6: 20BH
    if not (
        all(int(n) < int(nums[0]) and int(n) < int(nums[1]) for n in ['2', '0']) and  # both numbers too small
        'B' not in letters and 'H' not in letters and  # both letters wrong
        is_too_early('B', letters) and is_too_early('H', letters)  # both letters too early
    ):
        return "Failed 20BH"
    
    # Guess 7: 17XY
    if not (
        '7' in nums and nums[0] != '7' and  # 7 in wrong position
        'X' not in letters and 'Y' not in letters and  # both letters wrong
        is_too_late('X', letters) and is_too_late('Y', letters)  # both letters too late
    ):
        return "Failed 17XY"
    
    return "All passed"

# Test the combination
nums = ['1', '7']
letters = ['K', 'R']

result = verify_all_conditions(nums, letters)
print(result)