# Define the possible numbers and letters based on the deductions
possible_numbers = [6, 3]
possible_letters = ['O', 'Q']

# Check each condition to ensure the combination satisfies all
def check_combination(numbers, letters):
    # Convert numbers and letters to a string for easy comparison
    combination = f"{numbers[0]}{numbers[1]}{letters[0]}{letters[1]}"
    
    # Define the guesses and feedback
    guesses = [
        ("38UW", (1, 1, 0, 0)),  # one number correct but wrong position, one number incorrect and too large, both letters incorrect and too late
        ("69OG", (1, 1, 1, 1)),  # one number correct and correct position, one number incorrect and too large, one letter correct and correct position, one letter incorrect and too early
        ("83OP", (1, 1, 1, 0)),  # one number correct and correct position, one number incorrect and too large, one letter correct and correct position, one letter incorrect
        ("48SP", (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect
        ("42ZE", (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect
        ("01WG", (0, 0, 0, 0)),  # both numbers incorrect and too small, both letters incorrect
        ("07JX", (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect
        ("31HK", (1, 1, 0, 0)),  # one number correct but wrong position, one number incorrect and too small, both letters incorrect and too early
        ("37VY", (1, 1, 0, 0)),  # one number correct but wrong position, one number incorrect and too large, both letters incorrect and too late
        ("62QH", (1, 1, 0, 0)),  # one number correct and correct position, one number incorrect and too small, both letters incorrect
        ("94CL", (0, 0, 0, 0)),  # both numbers incorrect, both letters incorrect and too early
        ("71TE", (0, 0, 1, 0))   # both numbers incorrect, one letter correct but wrong position, one letter incorrect and too early
    ]
    
    # Check each guess against the combination
    for guess, feedback in guesses:
        num_correct_pos = sum(1 for i in range(2) if guess[i] == combination[i])
        num_correct_wrong_pos = sum(1 for i in range(2) if guess[i] in combination[:2] and guess[i] != combination[i])
        letter_correct_pos = sum(1 for i in range(2, 4) if guess[i] == combination[i])
        letter_correct_wrong_pos = sum(1 for i in range(2, 4) if guess[i] in combination[2:] and guess[i] != combination[i])
        
        if (num_correct_pos, num_correct_wrong_pos, letter_correct_pos, letter_correct_wrong_pos) != feedback:
            return False
    
    return True

# Check the combination
if check_combination(possible_numbers, possible_letters):
    print(f"The correct combination is: {possible_numbers[0]}{possible_numbers[1]}{possible_letters[0]}{possible_letters[1]}")
else:
    print("The combination does not satisfy all conditions.")