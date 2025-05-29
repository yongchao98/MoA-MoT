def check_guess(guess, actual):
    feedback = []
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            feedback.append("correct position")
        elif guess[i] in actual[:2]:
            feedback.append("correct wrong position")
        elif int(guess[i]) > int(actual[i]):
            feedback.append("too large")
        else:
            feedback.append("too small")
    
    # Check letters
    for i in range(2):
        if guess[i+2] == actual[i+2]:
            feedback.append("correct position")
        elif guess[i+2] in actual[2:]:
            feedback.append("correct wrong position")
        elif ord(guess[i+2]) > ord(actual[i+2]):
            feedback.append("too late")
        else:
            feedback.append("too early")
    return feedback

def matches_feedback(guess, actual, expected_feedback):
    feedback = check_guess(guess, actual)
    return feedback == expected_feedback

def find_password():
    numbers = [str(i).zfill(1) for i in range(10)]
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    # Store all guesses and their feedback
    guesses = {
        "87CF": ["too large", "too large", "too early", "too early"],
        "71SM": ["correct wrong position", "too large", "correct wrong position", "too late"],
        "23AY": ["incorrect", "incorrect", "incorrect", "incorrect"],
        "53PD": ["incorrect", "incorrect", "correct wrong position", "too early"],
        "01UF": ["correct wrong position", "too small", "incorrect", "incorrect"],
        "43US": ["incorrect", "incorrect", "too late", "too late"]
    }
    
    # Try all possible combinations
    for n1 in numbers:
        for n2 in numbers:
            if n1 != n2:  # Numbers can't repeat
                for l1 in letters:
                    for l2 in letters:
                        if l1 != l2:  # Letters can't repeat
                            candidate = n1 + n2 + l1 + l2
                            valid = True
                            
                            # Check if candidate matches all given feedback
                            for guess, feedback in guesses.items():
                                if not matches_feedback(guess, candidate, feedback):
                                    valid = False
                                    break
                            
                            if valid:
                                return [n1, n2, l1, l2]
    return None

result = find_password()
print(result)