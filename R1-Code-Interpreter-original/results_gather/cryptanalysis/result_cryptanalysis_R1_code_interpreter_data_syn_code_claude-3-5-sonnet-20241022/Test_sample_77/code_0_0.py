def check_guess(guess, actual):
    # Returns feedback similar to what we have in the conditions
    feedback = []
    for i in range(2):  # Check numbers
        if guess[i] == actual[i]:
            feedback.append("correct position")
        elif guess[i] in actual[:2]:
            feedback.append("correct wrong position")
        elif int(guess[i]) > int(actual[i]):
            feedback.append("too large")
        else:
            feedback.append("incorrect")
    
    for i in range(2):  # Check letters
        if guess[i+2] == actual[i+2]:
            feedback.append("correct position")
        elif guess[i+2] in actual[2:]:
            feedback.append("correct wrong position")
        else:
            feedback.append("incorrect")
    return feedback

def matches_condition(candidate):
    # Check against all known conditions
    conditions = [
        ("65WE", ["too large", "too large", "incorrect", "incorrect"]),
        ("82RV", ["correct wrong position", "too large", "incorrect", "incorrect"]),
        ("12KX", ["correct wrong position", "incorrect", "incorrect", "incorrect"]),
        ("97FO", ["too large", "too large", "correct wrong position", "incorrect"]),
        ("38AC", ["too large", "too large", "incorrect", "incorrect"]),
        ("03CX", ["correct wrong position", "too large", "incorrect", "incorrect"]),
        ("85XY", ["too large", "too large", "correct wrong position", "incorrect"])
    ]
    
    for guess, expected_feedback in conditions:
        feedback = check_guess(guess, candidate)
        if feedback != expected_feedback:
            return False
    return True

# Generate all possible combinations
possible_numbers = [str(i).zfill(2) for i in range(100)]
possible_letters = [chr(i) for i in range(65, 91)]  # A-Z

found = None
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        numbers = f"{n1}{n2}"
        for l1 in possible_letters:
            for l2 in possible_letters:
                if l1 == l2:
                    continue
                candidate = numbers + l1 + l2
                if matches_condition(candidate):
                    found = candidate
                    break
            if found:
                break
        if found:
            break
    if found:
        break

if found:
    print([found[0], found[1], found[2], found[3]])