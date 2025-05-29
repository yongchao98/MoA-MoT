def check_guess(guess):
    # All conditions from the hints
    if guess[0] == '0' or guess[1] == '5': return False  # Hint 1
    if guess[0] == '3' and guess[1] == '9': return False  # Hint 6
    if guess[0] == '5' or guess[1] == '4': return False  # Hint 10
    
    # Check for letter L in correct position (third position)
    if guess[2] != 'L': return False
    
    # Check for invalid letters
    invalid_letters = set('XUCHPJOTKREABZV')
    if guess[2] in invalid_letters or guess[3] in invalid_letters:
        return False
    
    # Numbers must be different
    if guess[0] == guess[1]: return False
    
    # Letters must be different
    if guess[2] == guess[3]: return False
    
    return True

def find_password():
    numbers = '123456789'
    letters = 'FGHIJKLMNOPQRSTUVWXYZ'  # Starting after B (too early in alphabet)
    
    valid_combinations = []
    
    for n1 in numbers:
        for n2 in numbers:
            for l2 in letters:  # l1 must be 'L'
                guess = n1 + n2 + 'L' + l2
                if check_guess(guess):
                    valid_combinations.append([n1, n2, 'L', l2])
    
    return valid_combinations

result = find_password()
print(result)