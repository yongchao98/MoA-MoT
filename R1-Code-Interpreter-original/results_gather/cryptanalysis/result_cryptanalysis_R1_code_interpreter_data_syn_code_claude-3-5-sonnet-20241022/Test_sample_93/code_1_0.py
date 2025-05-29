def check_guess(guess, actual):
    # Returns feedback based on the guess
    correct_num = 0
    correct_pos_num = 0
    correct_let = 0
    correct_pos_let = 0
    
    # Check numbers
    for i in range(2):
        if guess[i] == actual[i]:
            correct_pos_num += 1
        elif guess[i] in actual[:2]:
            correct_num += 1
            
    # Check letters
    for i in range(2,4):
        if guess[i] == actual[i]:
            correct_pos_let += 1
        elif guess[i] in actual[2:]:
            correct_let += 1
    
    return (correct_pos_num, correct_num, correct_pos_let, correct_let)

def is_valid_password(password):
    # Check against all known guesses
    guesses = [
        ("15TV", (0,0,0,0)),  # both numbers wrong, both letters too late
        ("07LC", (0,0,0,0)),  # both numbers wrong, both letters wrong
        ("63ZT", (0,0,0,0)),  # both numbers wrong, both letters too late
        ("02MJ", (0,0,0,0)),  # both numbers too small, both letters wrong
        ("08RE", (0,0,0,0)),  # both numbers wrong, both letters wrong
        ("23NC", (0,0,0,0)),  # both numbers too small, both letters wrong
        ("54GB", (1,0,0,0)),  # one number correct position, one wrong
        ("20DI", (0,0,0,1)),  # numbers too small, one letter wrong position
        ("74VO", (1,0,0,1))   # one number correct position, one letter wrong position
    ]
    
    for guess, feedback in guesses:
        result = check_guess(guess, password)
        if result != feedback:
            return False
            
    return True

# Generate all possible combinations
possible_numbers = [4,5,7,9]  # Numbers that could be in the solution
possible_letters = [chr(x) for x in range(65,91) if x not in [67,68,69,76,82,84,86,90]]  # Excluding C,D,E,L,R,T,V,Z

for n1 in possible_numbers:
    for n2 in possible_numbers:
        if n1 != n2:  # Numbers can't repeat
            for l1 in possible_letters:
                for l2 in possible_letters:
                    if l1 != l2:  # Letters can't repeat
                        password = f"{n1}{n2}{l1}{l2}"
                        if is_valid_password(password):
                            print(f"['{n1}', '{n2}', '{l1}', '{l2}']")