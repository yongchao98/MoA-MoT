from itertools import permutations
import string

def check_numbers(guess_nums, actual_nums, feedback):
    if feedback == "both incorrect":
        return all(n not in actual_nums for n in guess_nums)
    elif feedback == "both incorrect and too small":
        return all(n not in actual_nums and int(n) < min(int(x) for x in actual_nums) for n in guess_nums)
    elif feedback == "one correct in position":
        return sum(g == a for g, a in zip(guess_nums, actual_nums)) == 1
    elif feedback == "one correct wrong position":
        return sum(g in actual_nums for g in guess_nums) == 1 and all(g != a for g, a in zip(guess_nums, actual_nums))
    return False

def check_letters(guess_letters, actual_letters, feedback):
    if feedback == "both incorrect":
        return all(l not in actual_letters for l in guess_letters)
    elif feedback == "one correct wrong position":
        return sum(l in actual_letters for l in guess_letters) == 1 and all(g != a for g, a in zip(guess_letters, actual_letters))
    return False

def is_valid_combination(nums, letters):
    # Check guess 1: 84KP
    if not (check_numbers(['8', '4'], nums, "both incorrect") and
            check_letters(['K', 'P'], letters, "one correct wrong position")):
        return False
    
    # Check guess 2: 13GY
    if not (check_numbers(['1', '3'], nums, "one correct in position") and
            check_letters(['G', 'Y'], letters, "both incorrect")):
        return False
    
    # Check guess 3: 31HM
    if not (check_numbers(['3', '1'], nums, "one correct wrong position") and
            check_letters(['H', 'M'], letters, "both incorrect")):
        return False
    
    # Check guess 4: 84RL
    if not (check_numbers(['8', '4'], nums, "both incorrect") and
            check_letters(['R', 'L'], letters, "one correct wrong position")):
        return False
    
    # Check guess 5: 02AN
    if not (check_numbers(['0', '2'], nums, "both incorrect and too small") and
            check_letters(['A', 'N'], letters, "both incorrect")):
        return False
    
    # Check guess 6: 20BH
    if not (check_numbers(['2', '0'], nums, "both incorrect and too small") and
            check_letters(['B', 'H'], letters, "both incorrect")):
        return False
    
    # Check guess 7: 17XY
    if not (check_numbers(['1', '7'], nums, "one correct wrong position") and
            check_letters(['X', 'Y'], letters, "both incorrect")):
        return False
    
    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = list(string.ascii_uppercase)

for nums in permutations(numbers, 2):
    for lets in permutations(letters, 2):
        if is_valid_combination(nums, lets):
            print([nums[0], nums[1], lets[0], lets[1]])