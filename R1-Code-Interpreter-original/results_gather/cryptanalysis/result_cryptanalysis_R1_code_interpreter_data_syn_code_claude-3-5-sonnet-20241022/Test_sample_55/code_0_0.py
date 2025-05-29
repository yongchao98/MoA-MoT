from itertools import permutations
import string

def is_consistent_with_feedback(password, guesses_feedback):
    # Convert password to list for easier comparison
    password = list(password)
    
    for guess, feedback in guesses_feedback:
        guess = list(guess)
        # Check each feedback condition
        if guess == password:
            return False  # We know none of the guesses is correct
            
        # Count correct numbers and letters
        correct_nums = 0
        correct_letters = 0
        wrong_pos_nums = 0
        wrong_pos_letters = 0
        
        # First check for exact matches
        for i in range(4):
            if i < 2:  # Numbers
                if guess[i] == password[i]:
                    correct_nums += 1
            else:  # Letters
                if guess[i] == password[i]:
                    correct_letters += 1
        
        # Then check for wrong positions
        for i in range(2):
            if guess[i] in password[:2] and guess[i] != password[i]:
                wrong_pos_nums += 1
        for i in range(2,4):
            if guess[i] in password[2:] and guess[i] != password[i]:
                wrong_pos_letters += 1
        
        # Verify against specific feedback
        if guess == ['4','3','A','J']:
            if wrong_pos_nums != 1 or correct_nums != 0:
                return False
            if not (ord(guess[2]) < ord(password[2]) and ord(guess[3]) < ord(password[3])):
                return False
        elif guess == ['4','0','H','Z']:
            if wrong_pos_nums != 1 or correct_nums != 0:
                return False
        elif guess == ['8','2','Y','N']:
            if correct_nums != 0 or correct_letters != 1:
                return False
            if not any(ord(guess[i]) < ord(password[i]) for i in [2,3]):
                return False
        elif guess == ['1','5','N','V']:
            if correct_nums != 0 or correct_letters != 0:
                return False
        elif guess == ['3','9','U','B']:
            if correct_nums != 0 or correct_letters != 0:
                return False
        elif guess == ['1','8','U','J']:
            if correct_nums != 0 or correct_letters != 0:
                return False
        elif guess == ['2','9','F','P']:
            if correct_nums != 0:
                return False
            if not (ord(guess[2]) < ord(password[2]) and ord(guess[3]) < ord(password[3])):
                return False
        elif guess == ['2','6','Q','L']:
            if correct_nums != 0 or (correct_letters + wrong_pos_letters) != 1:
                return False
            if not any(ord(guess[i]) < ord(password[i]) for i in [2,3]):
                return False
    
    return True

# Generate all possible combinations
numbers = '0123456789'
letters = string.ascii_uppercase

guesses = [
    ('43AJ', 'one number wrong pos, one too small, both letters too early'),
    ('40HZ', 'one number wrong pos, one too small'),
    ('82YN', 'both numbers wrong, one letter correct pos'),
    ('15NV', 'all wrong'),
    ('39UB', 'all wrong'),
    ('18UJ', 'all wrong'),
    ('29FP', 'both numbers wrong, both letters too early'),
    ('26QL', 'both numbers wrong, one letter wrong pos, one too early')
]

# Try all possible combinations
valid_passwords = []
for nums in permutations(numbers, 2):
    for lets in permutations(letters, 2):
        password = nums + lets
        if is_consistent_with_feedback(password, guesses):
            valid_passwords.append(list(password))

print(valid_passwords)