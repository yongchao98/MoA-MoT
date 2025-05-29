def analyze_guesses():
    # Store all guesses and their feedback
    guesses = [
        ("20VX", "both nums wrong & small, both letters wrong"),
        ("73GA", "both nums wrong, letters wrong & early"),
        ("13KM", "both nums wrong & small, letters wrong & early"),
        ("74CD", "one num correct pos, one wrong, letters wrong & early"),
        ("57GT", "both nums wrong, letters wrong & early"),
        ("32XO", "both nums wrong & small, both letters wrong"),
        ("82EL", "one num correct pos, one small, letters wrong & early"),
        ("25DA", "both nums wrong, letters wrong & early"),
        ("87WU", "one num correct pos, one wrong, one letter correct pos, one wrong"),
        ("95ZM", "both nums wrong, one letter correct pos, one early")
    ]
    
    # Analyze numbers
    possible_nums = set(range(0, 10))
    for guess, feedback in guesses:
        num1, num2 = int(guess[0]), int(guess[1])
        
        if "both nums wrong & small" in feedback:
            possible_nums -= set(range(0, max(num1, num2) + 1))
        elif "one num correct pos" in feedback:
            print(f"Guess with one correct number position: {guess}")
            
    print("\nPossible numbers remaining:", sorted(list(possible_nums)))
    
    # Analyze letters
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for guess, feedback in guesses:
        letters = guess[2:]
        if "one letter correct pos" in feedback:
            print(f"Guess with one correct letter position: {guess}")
            
analyze_guesses()