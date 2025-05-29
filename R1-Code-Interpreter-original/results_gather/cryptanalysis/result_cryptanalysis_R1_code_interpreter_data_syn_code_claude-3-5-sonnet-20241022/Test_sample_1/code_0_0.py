def check_guess(guess, actual):
    feedback = []
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            feedback.append("correct position")
        elif guess[i] in actual[:2]:
            feedback.append("wrong position")
        elif int(guess[i]) < int(actual[i]):
            feedback.append("too small")
        else:
            feedback.append("incorrect")
    
    # Check letters
    for i in range(2,4):
        if guess[i] == actual[i]:
            feedback.append("correct position")
        elif guess[i] in actual[2:]:
            feedback.append("wrong position")
        elif ord(guess[i]) > ord(actual[i]):
            feedback.append("too late")
        else:
            feedback.append("incorrect")
    return feedback

def matches_conditions(password):
    # Convert password to list of strings
    pwd = [str(password[0]), str(password[1]), password[2], password[3]]
    
    # Check all conditions
    g1 = check_guess(['9','2','X','L'], pwd)
    if not (g1.count("wrong position") == 1 and "too small" in g1 and g1.count("incorrect") == 2):
        return False
        
    g2 = check_guess(['7','5','T','E'], pwd)
    if not (g2.count("incorrect") == 4):
        return False
        
    g3 = check_guess(['3','8','Y','F'], pwd)
    if not (g3.count("correct position") == 1 and g3.count("incorrect") == 3):
        return False
        
    g4 = check_guess(['1','4','M','X'], pwd)
    if not (g4.count("incorrect") == 4):
        return False
        
    g5 = check_guess(['5','6','C','Y'], pwd)
    if not (g5.count("correct position") == 1 and "too late" in g5 and g5.count("incorrect") == 2):
        return False
        
    g6 = check_guess(['0','9','W','V'], pwd)
    if not (g6.count("correct position") == 1 and "too small" in g6 and 
            g6.count("wrong position") == 1 and g6.count("incorrect") == 1):
        return False
        
    return True

# Try all possible combinations
numbers = range(10)
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
solution = None

for n1 in numbers:
    for n2 in numbers:
        if n1 != n2:  # numbers can't repeat
            for l1 in letters:
                for l2 in letters:
                    if l1 != l2:  # letters can't repeat
                        if matches_conditions([n1, n2, l1, l2]):
                            solution = [n1, n2, l1, l2]
                            print(solution)