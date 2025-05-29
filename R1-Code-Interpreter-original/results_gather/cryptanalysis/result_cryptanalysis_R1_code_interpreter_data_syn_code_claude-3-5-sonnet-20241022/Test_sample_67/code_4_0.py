def verify_all_conditions(combo):
    def check_letters_position(guess_letter, feedback, combo_letter):
        if feedback == "too early":
            return ord(guess_letter) < ord(combo_letter)
        elif feedback == "too late":
            return ord(guess_letter) > ord(combo_letter)
        return True

    conditions = [
        # Format: [guess, number_feedback, letter1_feedback, letter2_feedback]
        ([7,9,'F','V'], "both too large", "incorrect", "incorrect"),
        ([3,2,'P','Z'], "both too small", "too late", "too late"),
        ([0,9,'E','F'], "both incorrect", "too early", "too early"),
        ([5,8,'Q','D'], "one wrong pos one large", "incorrect", "incorrect"),
        ([7,9,'O','J'], "both too large", "incorrect", "incorrect"),
        ([6,4,'L','Y'], "one wrong pos one large", "correct", "too late"),
        ([4,8,'H','I'], "one correct one large", "too early", "too early"),
        ([4,5,'T','G'], "both correct", "incorrect", "incorrect"),
        ([3,1,'I','B'], "both too small", "too early", "too early"),
        ([9,4,'V','W'], "one wrong pos one large", "too late", "too late"),
        ([7,0,'X','N'], "both incorrect", "too late", "too late"),
        ([7,0,'B','I'], "both incorrect", "too early", "too early"),
        ([8,9,'U','G'], "both too large", "incorrect", "incorrect"),
        ([7,0,'K','G'], "both incorrect", "too early", "too early")
    ]

    for condition in conditions:
        guess, num_feedback, l1_feedback, l2_feedback = condition
        
        # Check numbers
        if num_feedback == "both too large":
            if not (guess[0] > combo[0] and guess[1] > combo[1]):
                return False
        elif num_feedback == "both too small":
            if not (guess[0] < combo[0] and guess[1] < combo[1]):
                return False
        elif num_feedback == "both incorrect":
            if guess[0] in combo[:2] or guess[1] in combo[:2]:
                return False
        elif num_feedback == "both correct":
            if not (guess[0] == combo[0] and guess[1] == combo[1]):
                return False
        elif num_feedback == "one wrong pos one large":
            if not ((guess[0] in combo[:2] and guess[0] != combo[0] and guess[1] > max(combo[:2])) or
                    (guess[1] in combo[:2] and guess[1] != combo[1] and guess[0] > max(combo[:2]))):
                return False
        elif num_feedback == "one correct one large":
            if not (guess[0] == combo[0] and guess[1] > combo[1]):
                return False

        # Check letters
        if l1_feedback == "correct":
            if guess[2] != combo[2]:
                return False
        elif l1_feedback in ["too early", "too late"]:
            if not check_letters_position(guess[2], l1_feedback, combo[2]):
                return False
        elif l1_feedback == "incorrect":
            if guess[2] == combo[2] or guess[2] == combo[3]:
                return False

        if l2_feedback == "correct":
            if guess[3] != combo[3]:
                return False
        elif l2_feedback in ["too early", "too late"]:
            if not check_letters_position(guess[3], l2_feedback, combo[3]):
                return False
        elif l2_feedback == "incorrect":
            if guess[3] == combo[2] or guess[3] == combo[3]:
                return False

    return True

# Test all possible combinations
valid_solutions = []
for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    if letter != 'L':  # L is confirmed in third position
        candidate = [4, 5, 'L', letter]
        if verify_all_conditions(candidate):
            valid_solutions.append(candidate)

print("Solutions that satisfy ALL conditions:")
for solution in valid_solutions:
    print(solution)