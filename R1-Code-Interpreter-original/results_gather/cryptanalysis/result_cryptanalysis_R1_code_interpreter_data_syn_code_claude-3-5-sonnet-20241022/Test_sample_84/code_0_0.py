def analyze_guesses():
    # Initialize possible digits and letters
    possible_digits = set(range(10))
    possible_letters = set('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    
    # Known positions and exclusions
    first_digit_known = set()
    second_digit_known = set()
    first_letter_known = set()
    second_letter_known = set()
    
    # Analyze each guess
    # Guess 3: 43AG - 4 is correct in position 1
    first_digit_known.add(4)
    
    # Guess 12: 83LK - 8 is correct in position 1, K is in wrong position
    # This means 4 must be in position 1 since we can't have both 4 and 8 correct
    # K must be in first letter position
    first_letter_known.add('K')
    
    # Guess 17: 28IK - K is in wrong position
    # This confirms K must be in first letter position
    
    # From guesses with "too late in alphabet" feedback
    too_late = set('SWRQYUVWO')
    possible_letters -= too_late
    
    # From guesses with "too early in alphabet" feedback
    too_early = set('IDGE')
    possible_letters -= too_early
    
    # Numbers analysis
    # When 6 appears, it's in wrong position but correct
    # When 2 appears with feedback "one correct in wrong position"
    # This means 6 is the second digit
    second_digit_known.add(6)
    
    # Print results
    print("First digit must be:", 4)
    print("Second digit must be:", 6)
    print("First letter must be:", 'K')
    print("Possible second letters:", sorted(list(possible_letters - too_late - too_early - {'K'})))

analyze_guesses()