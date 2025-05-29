import string

def analyze_guesses():
    # All possible digits and letters
    digits = list(range(10))
    letters = list(string.ascii_uppercase)
    
    # From condition 8, we know the numbers are 4 and 5
    numbers = [4, 5]
    
    # Initialize possible letters
    possible_letters = set(letters)
    
    # From conditions about letters being too early or too late
    # Let's track what we know about letter positions
    for letter in 'FVPZEFQOJLYHITGIBVWXNBIKG':
        possible_letters.discard(letter)
    
    # From condition 6, one letter is correct and in position
    # Testing each position with remaining letters
    final_combinations = []
    
    # Try each remaining letter in each position
    for l1 in possible_letters:
        for l2 in possible_letters:
            if l1 != l2:  # Letters can't repeat
                # Test against condition 6 (64LY)
                if (l1 == 'L' or l2 == 'L'):  # One letter must be L
                    combination = [4, 5, l1, l2]
                    final_combinations.append(combination)
    
    print("Possible combinations:")
    for combo in final_combinations:
        print(combo)

analyze_guesses()