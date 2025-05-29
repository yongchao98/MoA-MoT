def analyze_first_digit():
    all_digits = set(range(10))
    impossible = set()
    
    # Add numbers we know are wrong
    for guess, feedback in [
        ('58TD', 'both incorrect'),
        ('19UI', 'both incorrect'),
        ('75AC', 'both incorrect'),
        ('51PB', 'both incorrect'),
        ('08LP', 'both incorrect')
    ]:
        impossible.add(int(guess[0]))
        impossible.add(int(guess[1]))
    
    # Add numbers that are too large
    for n in range(9, -1, -1):
        if n >= 9:  # from 98OJ feedback
            impossible.add(n)
            
    # Add numbers that are too small
    for n in range(0, 10):
        if n <= 0:  # from 02HU feedback
            impossible.add(n)
    
    possible = all_digits - impossible - {2}  # 2 must be in second position
    print("Possible first digits:", sorted(list(possible)))

analyze_first_digit()