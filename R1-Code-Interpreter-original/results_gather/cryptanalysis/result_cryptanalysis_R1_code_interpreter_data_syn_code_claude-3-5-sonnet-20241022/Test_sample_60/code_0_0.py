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
    for i in range(2,4):
        if guess[i] == actual[i]:
            feedback.append("correct position")
        elif guess[i] in actual[2:]:
            feedback.append("correct wrong position")
        elif ord(guess[i]) > ord(actual[i]):
            feedback.append("too late")
        else:
            feedback.append("too early")
    return feedback

def matches_conditions(password):
    # Convert password list to proper format for checking
    password_str = ''.join(str(x) for x in password)
    
    # Check all conditions
    g1 = check_guess("15PU", password_str)
    if not (g1.count("correct wrong position") == 1 and 
            g1.count("too large") == 1 and 
            g1.count("too late") == 2):
        return False

    g2 = check_guess("43EH", password_str)
    if not (g2.count("too large") == 2 and 
            g2.count("correct wrong position") == 1 and 
            g2.count("too early") == 1):
        return False

    g3 = check_guess("94HG", password_str)
    if not (g3.count("too large") == 2 and 
            g3.count("correct position") == 1 and 
            g3.count("too early") == 1):
        return False

    g4 = check_guess("36NW", password_str)
    if not (g4.count("too large") == 2 and 
            g4.count("too late") == 2):
        return False

    g5 = check_guess("39FU", password_str)
    if not (g5.count("too large") == 2 and 
            g5.count("incorrect") == 2):
        return False

    g6 = check_guess("27JZ", password_str)
    if not (g6.count("too large") == 2 and 
            g6.count("correct wrong position") == 1 and 
            g6.count("too late") == 1):
        return False

    return True

# Try all possible combinations
numbers = range(0, 10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        password = [str(n1), str(n2), l1, l2]
                        if matches_conditions(password):
                            print(password)