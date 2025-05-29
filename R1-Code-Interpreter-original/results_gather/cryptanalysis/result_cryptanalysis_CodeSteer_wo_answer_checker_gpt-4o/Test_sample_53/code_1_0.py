def is_valid_combination(numbers, letters, feedback):
    for guess, fb in feedback:
        guess_numbers, guess_letters = guess[:2], guess[2:]
        
        # Check numbers
        if fb['numbers'] == 'both incorrect and too small':
            if not all(n > int(g) for n, g in zip(numbers, guess_numbers)):
                return False
        elif fb['numbers'] == 'one correct but in wrong position':
            if not any(n == int(g) for n, g in zip(numbers, guess_numbers)):
                return False
        elif fb['numbers'] == 'both incorrect':
            if any(n == int(g) for n, g in zip(numbers, guess_numbers)):
                return False
        
        # Check letters
        if fb['letters'] == 'one correct but in wrong position':
            if not any(l in guess_letters for l in letters):
                return False
        elif fb['letters'] == 'both incorrect and too early':
            if not all(l > g for l, g in zip(letters, guess_letters)):
                return False
        elif fb['letters'] == 'both incorrect':
            if any(l == g for l, g in zip(letters, guess_letters)):
                return False

    return True

def find_password():
    possible_numbers = [4, 5, 6, 8]
    possible_letters = ['S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    
    feedback = [
        ("30AU", {'numbers': 'both incorrect and too small', 'letters': 'one correct but in wrong position'}),
        ("60KI", {'numbers': 'one correct but in wrong position', 'letters': 'both incorrect and too early'}),
        ("23QB", {'numbers': 'both incorrect and too small', 'letters': 'both incorrect'}),
        ("03HC", {'numbers': 'both incorrect and too small', 'letters': 'both incorrect and too early'}),
        ("79MI", {'numbers': 'both incorrect', 'letters': 'both incorrect'})
    ]
    
    for n1 in possible_numbers:
        for n2 in possible_numbers:
            if n1 == n2:
                continue
            for l1 in possible_letters:
                for l2 in possible_letters:
                    if l1 == l2:
                        continue
                    if is_valid_combination([n1, n2], [l1, l2], feedback):
                        return [str(n1), str(n2), l1, l2]

password = find_password()
print(password)