from itertools import product

def check_guess(guess, numbers, letters):
    # Returns True if guess is consistent with all conditions
    
    def check_condition(test_guess, feedback):
        g_nums = test_guess[:2]
        g_lets = test_guess[2:]
        
        if feedback == "both numbers incorrect":
            if numbers[0] in g_nums or numbers[1] in g_nums:
                return False
        elif feedback == "one correct position":
            correct_count = sum(1 for i in range(2) if numbers[i] == g_nums[i])
            if correct_count != 1:
                return False
        
        if feedback == "both letters incorrect":
            if letters[0] in g_lets or letters[1] in g_lets:
                return False
        elif feedback == "one letter correct position":
            correct_count = sum(1 for i in range(2) if letters[i] == g_lets[i])
            if correct_count != 1:
                return False
        
        return True

    # Check all conditions
    if not check_condition("87VX", "both numbers incorrect, both letters incorrect late"):
        return False
    if not check_condition("26AP", "both numbers incorrect, both letters incorrect"):
        return False
    if not check_condition("98ZT", "one number correct position, both letters incorrect late"):
        return False
    if not check_condition("98HL", "one number correct position, one letter correct position"):
        return False
    if not check_condition("93BN", "one number correct position, both letters incorrect"):
        return False
    if not check_condition("37CG", "both numbers incorrect, both letters incorrect early"):
        return False
    if not check_condition("17KW", "both numbers incorrect, both letters incorrect late"):
        return False
    if not check_condition("92XC", "one number correct position, both letters incorrect"):
        return False
    if not check_condition("30NX", "both numbers incorrect small, both letters incorrect late"):
        return False
    if not check_condition("74DJ", "one number correct position, both letters incorrect"):
        return False
    
    return True

# Generate all possible combinations
digits = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for nums in product(digits, repeat=2):
    if nums[0] == nums[1]:  # Skip repeating numbers
        continue
    for lets in product(letters, repeat=2):
        if lets[0] == lets[1]:  # Skip repeating letters
            continue
        current_guess = [nums[0], nums[1], lets[0], lets[1]]
        if check_guess(current_guess, nums, lets):
            print(current_guess)