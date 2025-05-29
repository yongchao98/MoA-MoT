def check_number_relation(num, target):
    return "correct" if num == target else "too_small" if num < target else "too_large"

def check_letter_relation(letter, target):
    return "correct" if letter == target else "too_early" if letter < target else "too_late"

def validate_guess(guess, actual):
    # Compare numbers
    num_feedback = []
    for i in range(2):
        if guess[i] == actual[i]:
            num_feedback.append(("correct", "correct_pos"))
        elif guess[i] == actual[1-i]:
            num_feedback.append(("correct", "wrong_pos"))
        else:
            num_feedback.append((check_number_relation(int(guess[i]), int(actual[i])), "incorrect"))
    
    # Compare letters
    letter_feedback = []
    for i in range(2,4):
        if guess[i] == actual[i]:
            letter_feedback.append(("correct", "correct_pos"))
        elif guess[i] == actual[5-i]:
            letter_feedback.append(("correct", "wrong_pos"))
        else:
            letter_feedback.append((check_letter_relation(guess[i], actual[i]), "incorrect"))
    
    return num_feedback, letter_feedback

def matches_condition(candidate, guess, feedback):
    test_num_fb, test_letter_fb = validate_guess(guess, candidate)
    
    # Convert feedback to comparable format
    if feedback[0] == "both incorrect":
        expected_num = [("incorrect", "incorrect"), ("incorrect", "incorrect")]
    elif feedback[0] == "one correct wrong pos":
        expected_num = [("correct", "wrong_pos"), ("incorrect", None)]
    elif feedback[0] == "one correct correct pos":
        expected_num = [("correct", "correct_pos"), ("incorrect", None)]
    
    if feedback[1] == "both incorrect":
        expected_letter = [("incorrect", "incorrect"), ("incorrect", "incorrect")]
    elif feedback[1] == "one correct correct pos":
        expected_letter = [("correct", "correct_pos"), ("incorrect", None)]
    elif feedback[1] == "one correct wrong pos":
        expected_letter = [("correct", "wrong_pos"), ("incorrect", None)]
    
    # Additional conditions for too small/large
    if "too small" in feedback:
        for i, (status, pos) in enumerate(test_num_fb):
            if status != "correct" and not (status == "too_small") == ("too small" in feedback):
                return False
    
    if "too late" in feedback or "too early" in feedback:
        for i, (status, pos) in enumerate(test_letter_fb):
            if status != "correct":
                if "too late" in feedback and status != "too_late":
                    return False
                if "too early" in feedback and status != "too_early":
                    return False
    
    return True

def find_password():
    conditions = [
        ("65HO", ("one correct wrong pos", "both incorrect", "too small", None)),
        ("47IM", ("both incorrect", "both incorrect")),
        ("10GI", ("both incorrect", "both incorrect", "too small", "too early")),
        ("70XG", ("both incorrect", "both incorrect")),
        ("12JU", ("both incorrect", "one correct correct pos", "too small", "too late")),
        ("52SZ", ("both incorrect", "both incorrect", "too small", "too late")),
        ("64AV", ("one correct wrong pos", "both incorrect", "too small")),
        ("85BO", ("one correct correct pos", "both incorrect", "too small")),
        ("98RX", ("one correct wrong pos", "both incorrect", "too large", "too late")),
        ("60JU", ("one correct wrong pos", "one correct correct pos", "too small", "too late")),
        ("75UC", ("both incorrect", "both incorrect")),
        ("79NC", ("both incorrect", "one correct wrong pos", "too early"))
    ]
    
    # Try all possible combinations
    numbers = range(0, 100)
    letters = [chr(i) for i in range(65, 91)]  # A-Z
    
    for n1 in range(10):
        for n2 in range(10):
            if n1 == n2:
                continue
            for l1 in letters:
                for l2 in letters:
                    if l1 == l2:
                        continue
                    candidate = f"{n1}{n2}{l1}{l2}"
                    valid = True
                    for guess, feedback in conditions:
                        if not matches_condition(candidate, guess, feedback):
                            valid = False
                            break
                    if valid:
                        print(f"Found solution: {candidate}")

find_password()