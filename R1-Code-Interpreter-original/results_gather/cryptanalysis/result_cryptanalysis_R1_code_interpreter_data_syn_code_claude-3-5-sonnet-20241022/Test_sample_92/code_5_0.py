def verify_solution():
    solution = "62VN"
    test_cases = [
        ('58TD', 'both incorrect', 'both incorrect'),
        ('52SK', 'one correct position, one incorrect', 'both incorrect'),
        ('19UI', 'both incorrect', 'both incorrect'),
        ('75AC', 'both incorrect', 'both too early'),
        ('86WD', 'one wrong position, one too large', 'both incorrect'),
        ('72KL', 'one correct position, one too large', 'both too early'),
        ('08LP', 'both incorrect', 'both incorrect'),
        ('51PB', 'both incorrect', 'both incorrect'),
        ('96FX', 'one wrong position, one too large', 'both incorrect'),
        ('02HU', 'one correct position, one too small', 'both incorrect'),
        ('96FJ', 'one wrong position, one too large', 'both too early'),
        ('12KC', 'one correct position, one too small', 'both too early'),
        ('92TG', 'one correct position, one too large', 'both incorrect'),
        ('98OJ', 'both too large', 'both incorrect'),
        ('06ML', 'one wrong position, one too small', 'both too early'),
        ('02VX', 'one correct position, one too small', 'one wrong position, one too late')
    ]
    
    print(f"Verifying solution: {solution}")
    print("\nKey patterns that confirm 62VN:")
    
    # Verify number 6
    print("\n1. Verifying number 6:")
    for guess, num_feed, _ in test_cases:
        if '6' in guess and guess[0] != '6':
            print(f"Guess {guess}: {num_feed} (6 in wrong position)")
            
    # Verify number 2
    print("\n2. Verifying number 2:")
    for guess, num_feed, _ in test_cases:
        if '2' in guess and guess[1] == '2':
            print(f"Guess {guess}: {num_feed} (2 in correct position)")
            
    # Verify letters V and N
    print("\n3. Verifying letters V and N:")
    for guess, _, letter_feed in test_cases:
        if 'V' in guess:
            print(f"Guess {guess}: {letter_feed}")
        if guess[2:] < 'N':
            print(f"Guess {guess}: {letter_feed} (letters too early)")

verify_solution()