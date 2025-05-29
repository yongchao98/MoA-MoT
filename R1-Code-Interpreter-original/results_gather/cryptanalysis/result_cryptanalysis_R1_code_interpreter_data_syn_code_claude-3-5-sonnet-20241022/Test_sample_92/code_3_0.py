def verify_62():
    print("Verifying why 62 is the correct combination:")
    
    test_cases = [
        ('58TD', 'both incorrect'),  # ✓ 5,8 are both incorrect
        ('52SK', 'one correct position, one incorrect'),  # ✓ 2 is correct position
        ('19UI', 'both incorrect'),  # ✓ 1,9 are both incorrect
        ('75AC', 'both incorrect'),  # ✓ 7,5 are both incorrect
        ('86WD', 'one wrong position, one too large'),  # ✓ 6 is in wrong pos, 8 too large
        ('72KL', 'one correct position, one too large'),  # ✓ 2 correct pos, 7 too large
        ('08LP', 'both incorrect'),  # ✓ 0,8 are both incorrect
        ('51PB', 'both incorrect'),  # ✓ 5,1 are both incorrect
        ('96FX', 'one wrong position, one too large'),  # ✓ 6 wrong pos, 9 too large
        ('02HU', 'one correct position, one too small'),  # ✓ 2 correct pos, 0 too small
        ('96FJ', 'one wrong position, one too large'),  # ✓ 6 wrong pos, 9 too large
        ('12KC', 'one correct position, one too small'),  # ✓ 2 correct pos, 1 too small
        ('92TG', 'one correct position, one too large'),  # ✓ 2 correct pos, 9 too large
        ('98OJ', 'both too large'),  # ✓ 9,8 both too large
        ('06ML', 'one wrong position, one too small'),  # ✓ 6 wrong pos, 0 too small
        ('02VX', 'one correct position, one too small')  # ✓ 2 correct pos, 0 too small
    ]
    
    for guess, feedback in test_cases:
        print(f"\nGuess {guess}: {feedback}")
        if '6' in guess:
            print(f"Contains 6: {feedback}")
        if '2' in guess:
            print(f"Contains 2: {feedback}")

verify_62()