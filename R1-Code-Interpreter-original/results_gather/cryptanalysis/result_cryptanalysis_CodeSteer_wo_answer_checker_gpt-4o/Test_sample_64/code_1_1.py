def is_valid_combination(numbers, letters, combination):
    # Check each guess against the feedback
    guesses = [
        ("46DL", (0, 1, 1, 0)),
        ("31GM", (1, 0, 0, 1)),
        ("75WL", (0, 1, 0, 1)),
        ("74IA", (0, 0, 0, 0)),
        ("04OB", (0, 0, 0, 0)),
        ("38FP", (1, 1, 0, 0))
    ]
    
    for guess, feedback in guesses:
        num_correct = sum(1 for i in range(2) if combination[i] == guess[i])
        letter_correct = sum(1 for i in range(2, 4) if combination[i] == guess[i])
        
        if (num_correct, letter_correct) != feedback[:2]:
            return False
        
        if feedback[2] == 1 and combination[2] != guess[2]:
            return False
        
        if feedback[3] == 1 and combination[3] != guess[3]:
            return False
    
    return True

def find_password():
    numbers = [str(i) for i in range(10)]
    letters = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
    
    # Based on feedback, we can eliminate some numbers and letters early
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
                    if is_valid_combination(numbers, letters, combination):
                        return combination

password = find_password()
print(f"<<< {password} >>>")