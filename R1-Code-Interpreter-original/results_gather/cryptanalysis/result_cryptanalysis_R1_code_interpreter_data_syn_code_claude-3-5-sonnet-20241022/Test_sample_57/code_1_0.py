from itertools import permutations

def check_guess(guess, actual):
    # Check numbers
    num_feedback = []
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
    letter_feedback = []
    for i in range(2,4):
        if guess[i] == actual[i]:
            letter_feedback.append("correct position")
        elif guess[i] in actual[2:]:
            letter_feedback.append("wrong position")
        elif ord(guess[i]) > ord(actual[i]):
            letter_feedback.append("too late")
        else:
            letter_feedback.append("incorrect")
    
    return num_feedback, letter_feedback

def matches_conditions(password):
    # Convert password list to proper format
    password_str = password[0] + password[1] + password[2] + password[3]
    
    # Check all conditions
    # Condition 1: 43HI
    num_fb1, let_fb1 = check_guess("43HI", password_str)
    if not (("wrong position" in num_fb1 and "too small" in num_fb1) and 
            let_fb1.count("incorrect") == 2):
        return False

    # Condition 2: 45QH
    num_fb2, let_fb2 = check_guess("45QH", password_str)
    if not (("wrong position" in num_fb2 and "incorrect" in num_fb2) and 
            let_fb2.count("incorrect") == 2):
        return False

    # Condition 3: 01ZA
    num_fb3, let_fb3 = check_guess("01ZA", password_str)
    if not (num_fb3.count("too small") == 2 and 
            ("wrong position" in let_fb3) and ("too late" in let_fb3)):
        return False

    # Condition 4: 64RO
    num_fb4, let_fb4 = check_guess("64RO", password_str)
    if not (("correct position" in num_fb4 and "incorrect" in num_fb4) and 
            let_fb4.count("too late") == 2):
        return False

    # Condition 5: 92SB
    num_fb5, let_fb5 = check_guess("92SB", password_str)
    if not (("correct position" in num_fb5 and "too small" in num_fb5) and 
            let_fb5.count("incorrect") == 2):
        return False

    # Condition 6: 46DY
    num_fb6, let_fb6 = check_guess("46DY", password_str)
    if not (("wrong position" in num_fb6 and "incorrect" in num_fb6) and 
            let_fb6.count("incorrect") == 2):
        return False

    # Condition 7: 81CM
    num_fb7, let_fb7 = check_guess("81CM", password_str)
    if not (num_fb7.count("incorrect") == 2 and 
            ("correct position" in let_fb7) and ("incorrect" in let_fb7)):
        return False

    return True

# Generate all possible combinations
numbers = [str(i) for i in range(10)]
letters = [chr(i) for i in range(65, 91)]  # A-Z

# Try all possible combinations
for num1, num2 in permutations(numbers, 2):
    for let1, let2 in permutations(letters, 2):
        password = [num1, num2, let1, let2]
        if matches_conditions(password):
            print(password)