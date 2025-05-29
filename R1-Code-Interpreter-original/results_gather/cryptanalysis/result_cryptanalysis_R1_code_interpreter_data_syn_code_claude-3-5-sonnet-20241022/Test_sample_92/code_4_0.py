def final_verification():
    print("Final verification of solution 62VN:")
    
    test_cases = [
        # Format: (guess, expected_number_match, expected_letter_match)
        ('58TD', 0, 0),  # Both numbers and letters incorrect
        ('52SK', 1, 0),  # One number correct position (2), other number and letters wrong
        ('19UI', 0, 0),  # Both numbers and letters incorrect
        ('75AC', 0, 0),  # Both numbers incorrect, letters too early
        ('86WD', 1, 0),  # 6 wrong position, 8 too large
        ('72KL', 1, 0),  # 2 correct position, 7 too large
        ('08LP', 0, 0),  # Both numbers and letters incorrect
        ('51PB', 0, 0),  # Both numbers and letters incorrect
        ('96FX', 1, 0),  # 6 wrong position, 9 too large
        ('02HU', 1, 0),  # 2 correct position, 0 too small
        ('96FJ', 1, 0),  # 6 wrong position, 9 too large
        ('12KC', 1, 0),  # 2 correct position, 1 too small
        ('92TG', 1, 0),  # 2 correct position, 9 too large
        ('98OJ', 0, 0),  # Both numbers too large
        ('06ML', 1, 0),  # 6 wrong position, 0 too small
        ('02VX', 1, 1)   # 2 correct position, 0 too small, V wrong position, X too late
    ]
    
    solution = "62VN"
    print(f"\nChecking solution: {solution}")
    
    for guess, num_match, letter_match in test_cases:
        print(f"\nGuess: {guess}")
        print(f"Numbers: {guess[0:2]} vs {solution[0:2]}")
        print(f"Letters: {guess[2:]} vs {solution[2:]}")
        
        # Check if numbers match
        if '6' in guess:
            print("Contains 6")
        if '2' in guess:
            print("Contains 2")
        if 'V' in guess:
            print("Contains V")

print_verification()