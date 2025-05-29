def is_valid_combination(combination):
    guesses = [
        ("46DL", (0, 1, 1, 0)),
        ("31GM", (1, 0, 0, 1)),
        ("75WL", (0, 1, 1, 0)),
        ("74IA", (0, 0, 0, 0)),
        ("04OB", (0, 0, 0, 0)),
        ("38FP", (1, 1, 0, 0))
    ]
    
    for guess, feedback in guesses:
        num_correct = sum(1 for i in range(2) if combination[i] == guess[i])
        num_wrong_position = sum(1 for i in range(2) if combination[i] in guess[:2] and combination[i] != guess[i])
        
        letter_correct = sum(1 for i in range(2, 4) if combination[i] == guess[i])
        letter_wrong_position = sum(1 for i in range(2, 4) if combination[i] in guess[2:] and combination[i] != guess[i])
        
        if (num_correct, num_wrong_position) != feedback[:2]:
            return False
        
        if (letter_correct, letter_wrong_position) != feedback[2:]:
            return False
    
    return True

def find_password():
    possible_numbers = ['3', '8']
    possible_letters = ['L', 'M']
    
    for num1 in possible_numbers:
        for num2 in possible_numbers:
            if num1 == num2:
                continue
            for letter1 in possible_letters:
                for letter2 in possible_letters:
                    if letter1 == letter2:
                        continue
                    combination = [num1, num2, letter1, letter2]
                    if is_valid_combination(combination):
                        return combination

password = find_password()
print(f"<<< {password} >>>")