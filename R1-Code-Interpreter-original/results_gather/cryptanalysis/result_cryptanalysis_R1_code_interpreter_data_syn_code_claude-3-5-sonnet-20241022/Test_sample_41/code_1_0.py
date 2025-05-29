def verify_guess(guess, target):
    # Returns feedback based on the guess compared to target
    n1, n2, l1, l2 = guess
    tn1, tn2, tl1, tl2 = target
    
    # Convert numbers to integers for comparison
    n1, n2 = int(n1), int(n2)
    tn1, tn2 = int(tn1), int(tn2)
    
    feedback = []
    
    # Check numbers
    if n1 == tn1 and n2 == tn2:
        feedback.append("both numbers correct position")
    elif n1 == tn1:
        feedback.append("first number correct position")
    elif n2 == tn2:
        feedback.append("second number correct position")
    elif n1 == tn2 or n2 == tn1:
        feedback.append("one number correct wrong position")
    else:
        feedback.append("both numbers incorrect")
    
    # Check letters
    if l1 == tl1 and l2 == tl2:
        feedback.append("both letters correct position")
    elif l1 == tl1:
        feedback.append("first letter correct position")
    elif l2 == tl2:
        feedback.append("second letter correct position")
    elif l1 == tl2 or l2 == tl1:
        feedback.append("one letter correct wrong position")
    else:
        feedback.append("both letters incorrect")
    
    return feedback

candidates = [['7', '3', 'D', 'M'], ['9', '0', 'D', 'M']]

# Test each candidate against known guesses
for candidate in candidates:
    print(f"\nTesting candidate: {candidate}")
    
    # Test against guess 3: 93CH
    feedback = verify_guess(['9', '3', 'C', 'H'], candidate)
    print(f"Guess 93CH feedback: {feedback}")
    
    # Test against guess 5: 07DM
    feedback = verify_guess(['0', '7', 'D', 'M'], candidate)
    print(f"Guess 07DM feedback: {feedback}")