import itertools
import string

def check_guess(guess, actual):
    # Check numbers
    n1, n2 = guess[:2]
    a1, a2 = actual[:2]
    
    # Check letters
    l1, l2 = guess[2:]
    al1, al2 = actual[2:]
    
    return [n1 == a1, n2 == a2, l1 == al1, l2 == al2]

def is_valid_password(password):
    # Convert password to string format for comparison
    pwd = [str(x) for x in password]
    
    # Condition 1: 74SE
    c1 = sum(check_guess("74SE", pwd)) == 1  # one letter correct position
    
    # Condition 4: 28HI
    c4 = check_guess("28HI", pwd)[0]  # first number correct position
    
    # Condition 8: 12YB - both numbers too small
    c8 = int(pwd[0]) > 1 and int(pwd[1]) > 2
    
    # Condition 9: 87WE - one number correct wrong position
    c9 = (pwd[0] == '7' or pwd[1] == '8') and not (pwd[0] == '8' or pwd[1] == '7')
    
    # Condition 12: 30VF - both numbers too small
    c12 = int(pwd[0]) > 3 and int(pwd[1]) > 0
    
    # Condition 13: 68QG
    c13 = (check_guess("68QG", pwd)[0] or check_guess("68QG", pwd)[1])  # one number correct position
    
    return c1 and c4 and c8 and c9 and c12 and c13

# Generate possible combinations
numbers = range(0, 10)
letters = string.ascii_uppercase

valid_passwords = []
for n1, n2 in itertools.permutations(numbers, 2):
    for l1, l2 in itertools.permutations(letters, 2):
        password = [str(n1), str(n2), l1, l2]
        if is_valid_password(password):
            valid_passwords.append(password)

print(valid_passwords)