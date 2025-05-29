from itertools import permutations

def check_guess(guess, actual):
    # Returns feedback matching the format in the clues
    num_feedback = []
    letter_feedback = []
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            num_feedback.append("correct position")
        elif guess[i] in actual[:2]:
            num_feedback.append("wrong position")
        elif int(guess[i]) < int(actual[i]):
            num_feedback.append("too small")
        else:
            num_feedback.append("incorrect")
            
    # Check letters
    for i in range(2,4):
        if guess[i] == actual[i]:
            letter_feedback.append("correct position")
        elif guess[i] in actual[2:]:
            letter_feedback.append("wrong position")
        else:
            letter_feedback.append("incorrect")
            
    return num_feedback, letter_feedback

def matches_clue(candidate, guess, num_feedback, letter_feedback):
    test_num, test_letter = check_guess(guess, candidate)
    
    # Count feedback types
    num_correct_pos = sum(1 for x in test_num if x == "correct position")
    num_wrong_pos = sum(1 for x in test_num if x == "wrong position")
    num_too_small = sum(1 for x in test_num if x == "too small")
    
    letter_correct_pos = sum(1 for x in test_letter if x == "correct position")
    letter_wrong_pos = sum(1 for x in test_letter if x == "wrong position")
    
    # Compare with given feedback
    if num_feedback == "both incorrect and too small":
        return num_correct_pos == 0 and num_too_small == 2
    elif num_feedback == "both incorrect":
        return num_correct_pos == 0 and num_wrong_pos == 0
    elif num_feedback == "one correct but wrong position, one too small":
        return num_wrong_pos == 1 and num_too_small == 1
    elif num_feedback == "one correct position, one too small":
        return num_correct_pos == 1 and num_too_small == 1
        
    if letter_feedback == "one wrong position, one incorrect":
        return letter_correct_pos == 0 and letter_wrong_pos == 1
    elif letter_feedback == "both incorrect":
        return letter_correct_pos == 0 and letter_wrong_pos == 0
    elif letter_feedback == "one correct position, one incorrect":
        return letter_correct_pos == 1 and letter_wrong_pos == 0
    
    return False

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z
candidates = []

for n1, n2 in permutations(numbers, 2):
    for l1, l2 in permutations(letters, 2):
        candidate = [n1, n2, l1, l2]
        
        # Check against all clues
        if (matches_clue(candidate, ['5','4','V','U'], "both incorrect and too small", "one wrong position, one incorrect") and
            matches_clue(candidate, ['1','9','L','A'], "both incorrect", "one wrong position, one incorrect") and
            matches_clue(candidate, ['3','8','J','V'], "one correct but wrong position, one too small", "one correct position, one incorrect") and
            matches_clue(candidate, ['1','4','E','C'], "both incorrect and too small", "both incorrect") and
            matches_clue(candidate, ['4','7','Q','U'], "one correct position, one too small", "both incorrect")):
            candidates.append(candidate)

print(candidates)