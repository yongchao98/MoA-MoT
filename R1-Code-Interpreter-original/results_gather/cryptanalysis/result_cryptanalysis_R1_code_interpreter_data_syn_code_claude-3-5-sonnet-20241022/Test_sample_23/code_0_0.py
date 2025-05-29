def is_consistent(password, guesses_feedback):
    # Helper function to check if a letter is too late in alphabet
    def is_too_late(guess_letter, actual_letter):
        return ord(guess_letter) <= ord(actual_letter)
    
    # Helper function to check if a number is too large
    def is_too_large(guess_num, actual_num):
        return int(guess_num) <= int(actual_num)
    
    for guess, feedback in guesses_feedback:
        # Convert password to string for easier comparison
        pwd = ''.join(map(str, password))
        
        # Count correct positions and values
        correct_pos_nums = sum(1 for i in range(2) if pwd[i] == guess[i])
        correct_pos_letters = sum(1 for i in range(2,4) if pwd[i] == guess[i])
        
        # Count correct values in wrong positions
        correct_val_nums = sum(1 for i in range(2) for j in range(2) 
                             if i != j and pwd[i] == guess[j])
        correct_val_letters = sum(1 for i in range(2,4) for j in range(2,4) 
                                if i != j and pwd[i] == guess[j])
        
        # Check each guess against its feedback
        if guess == "98LZ" and feedback == (1,1,2,0):
            if not (correct_pos_nums == 1 and 
                   sum(not is_too_large(g, p) for g,p in zip(guess[:2], pwd[:2])) == 1 and
                   correct_pos_letters == 0 and
                   sum(not is_too_late(g, p) for g,p in zip(guess[2:], pwd[2:])) == 2):
                return False
        
        elif guess == "82EM" and feedback == (0,1,1,1):
            if not (correct_pos_nums == 0 and correct_val_nums == 1 and
                   correct_pos_letters == 0 and correct_val_letters == 1 and
                   sum(not is_too_late(g, p) for g,p in zip(guess[2:], pwd[2:])) == 1):
                return False
        
        elif guess == "36ZI" and feedback == (0,0,2,0):
            if not (correct_pos_nums == 0 and correct_val_nums == 0 and
                   correct_pos_letters == 0 and
                   sum(not is_too_late(g, p) for g,p in zip(guess[2:], pwd[2:])) == 2):
                return False
        
        elif guess == "21HR" and feedback == (0,1,1,1):
            if not (correct_pos_nums == 0 and correct_val_nums == 1 and
                   correct_pos_letters == 1 and
                   sum(not is_too_late(g, p) for g,p in zip(guess[2:], pwd[2:])) == 1):
                return False
    
    return True

# Generate all possible combinations
numbers = list(range(10))
letters = [chr(i) for i in range(65, 91)]  # A-Z

guesses_feedback = [
    ("98LZ", (1,1,2,0)),  # (correct_pos_nums, wrong_nums_too_large, wrong_letters_too_late, correct_pos_letters)
    ("82EM", (0,1,1,1)),  # (correct_pos_nums, correct_val_nums, correct_val_letters, wrong_letters_too_late)
    ("36ZI", (0,0,2,0)),  # (correct_pos_nums, correct_val_nums, wrong_letters_too_late, correct_pos_letters)
    ("21HR", (0,1,1,1))   # (correct_pos_nums, correct_val_nums, correct_pos_letters, wrong_letters_too_late)
]

# Find all valid combinations
for n1 in range(10):
    for n2 in range(10):
        if n1 == n2:
            continue
        for l1 in letters:
            for l2 in letters:
                if l1 == l2:
                    continue
                password = [str(n1), str(n2), l1, l2]
                if is_consistent(password, guesses_feedback):
                    print(password)