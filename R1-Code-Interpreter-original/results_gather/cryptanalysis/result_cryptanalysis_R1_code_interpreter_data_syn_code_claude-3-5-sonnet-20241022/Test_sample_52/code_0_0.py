from itertools import permutations

def is_consistent_with_feedback(password, guesses_feedback):
    for guess, feedback in guesses_feedback:
        # Convert password and guess to list for easier comparison
        pwd = list(password)
        gss = list(guess)
        
        # Check numbers (first two positions)
        correct_nums = 0
        misplaced_nums = 0
        for i in range(2):
            if pwd[i] == gss[i]:
                correct_nums += 1
            elif gss[i] in pwd[:2]:
                misplaced_nums += 1
                
        # Check letters (last two positions)
        correct_letters = 0
        misplaced_letters = 0
        for i in range(2, 4):
            if pwd[i] == gss[i]:
                correct_letters += 1
            elif gss[i] in pwd[2:]:
                misplaced_letters += 1
        
        # Verify against feedback
        if guess == "15NT":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "21OD":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "51BI":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "06GE":
            if correct_nums != 1 or correct_letters != 1:
                return False
        elif guess == "32FY":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "13UN":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "37DU":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "15IG":
            if correct_nums > 0 or correct_letters != 1 or misplaced_letters != 1:
                return False
        elif guess == "25IA":
            if correct_nums > 0 or correct_letters > 0:
                return False
        elif guess == "93PA":
            if misplaced_nums != 1 or correct_nums > 0 or correct_letters > 0:
                return False
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

guesses = [
    "15NT", "21OD", "51BI", "06GE", "32FY",
    "13UN", "37DU", "15IG", "25IA", "93PA"
]

# Test all possible combinations
valid_passwords = []
for num1 in numbers:
    for num2 in numbers:
        if num1 != num2:
            for let1 in letters:
                for let2 in letters:
                    if let1 != let2:
                        password = num1 + num2 + let1 + let2
                        if is_consistent_with_feedback(password, [(g, None) for g in guesses]):
                            valid_passwords.append(password)

print(valid_passwords)