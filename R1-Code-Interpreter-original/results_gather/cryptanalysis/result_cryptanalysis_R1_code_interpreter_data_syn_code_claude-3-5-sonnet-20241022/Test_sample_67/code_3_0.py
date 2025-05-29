def verify_combination(combo):
    def check_guess(guess, feedback):
        # Compare numbers
        if feedback.get("numbers") == "both incorrect and too large":
            if not (guess[0] > combo[0] and guess[1] > combo[1]):
                return False
        elif feedback.get("numbers") == "both incorrect and too small":
            if not (guess[0] < combo[0] and guess[1] < combo[1]):
                return False
        elif feedback.get("numbers") == "both incorrect":
            if guess[0] in combo[:2] or guess[1] in combo[:2]:
                return False
        elif feedback.get("numbers") == "both correct and in position":
            if not (guess[0] == combo[0] and guess[1] == combo[1]):
                return False
        elif feedback.get("numbers") == "one correct but wrong position, one too large":
            if not ((guess[0] in combo[:2] and guess[0] != combo[0] and guess[1] > max(combo[:2])) or
                    (guess[1] in combo[:2] and guess[1] != combo[1] and guess[0] > max(combo[:2]))):
                return False
        
        # Compare letters
        if feedback.get("letters") == "both incorrect":
            if guess[2] in combo[2:] or guess[3] in combo[2:]:
                return False
        elif feedback.get("letters") == "one correct in position, one incorrect too late":
            correct_found = False
            for i in [2, 3]:
                if guess[i] == combo[i]:
                    correct_found = True
            if not correct_found:
                return False
        
        return True

    conditions = [
        {"guess": [7,9,'F','V'], 
         "feedback": {"numbers": "both incorrect and too large", "letters": "both incorrect"}},
        {"guess": [3,2,'P','Z'],
         "feedback": {"numbers": "both incorrect and too small", "letters": "both incorrect and too late"}},
        {"guess": [0,9,'E','F'],
         "feedback": {"numbers": "both incorrect", "letters": "both incorrect and too early"}},
        {"guess": [5,8,'Q','D'],
         "feedback": {"numbers": "one correct but wrong position, one too large", "letters": "both incorrect"}},
        {"guess": [7,9,'O','J'],
         "feedback": {"numbers": "both incorrect and too large", "letters": "both incorrect"}},
        {"guess": [6,4,'L','Y'],
         "feedback": {"numbers": "one correct but wrong position, one too large", 
                     "letters": "one correct in position, one incorrect too late"}},
        {"guess": [4,8,'H','I'],
         "feedback": {"numbers": "one correct in position, one too large", "letters": "both incorrect and too early"}},
        {"guess": [4,5,'T','G'],
         "feedback": {"numbers": "both correct and in position", "letters": "both incorrect"}},
    ]

    for condition in conditions:
        if not check_guess(condition["guess"], condition["feedback"]):
            return False
    return True

# Test all possible combinations that match our basic criteria
solutions = []
for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    if letter != 'L':
        candidate = [4, 5, 'L', letter]
        if verify_combination(candidate):
            solutions.append(candidate)

print("Valid solutions after checking all conditions:")
for solution in solutions:
    print(solution)