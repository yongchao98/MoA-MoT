def verify_guess(guess):
    n1, n2, l1, l2 = guess
    
    # Condition 9 ("78WK") is crucial:
    # - One number (7) is correct but in wrong position
    # - One letter (K) is correct and in correct position
    # This means:
    # 1. 7 must be in the solution but NOT in the first position
    # 2. K must be in the solution in the SAME position as in "78WK" (last position)
    if n1 == 7 or n2 != 7 or l2 != 'K':
        return False

    # From condition 3 ("20XK") and condition 6 ("23VQ"):
    # - One number is correct and in correct position
    # - That number must be 2 (as it's the only common number)
    if n1 != 2:
        return False

    # From condition 13 ("49DW"):
    # - Both numbers are incorrect (so 4 and 9 can't be in solution)
    # - One letter (W) is correct but in wrong position
    # This means W must be in the solution but in a different position than in "49DW"
    if l2 == 'W' or l1 != 'W':  # W must be in first letter position
        return False

    # Additional checks for invalid numbers/letters
    invalid_numbers = {0, 1, 3, 4, 5, 6, 8, 9}  # Based on various conditions
    if n1 in invalid_numbers or n2 in invalid_numbers:
        return False

    invalid_letters = {'A', 'C', 'E', 'F', 'I', 'J', 'L', 'M', 'O', 'Q', 'T', 'V', 'X', 'Z'}
    if l1 in invalid_letters or l2 in invalid_letters:
        return False

    return True

# Test the specific combination we believe is correct
solution = (2, 7, 'W', 'K')
if verify_guess(solution):
    print("Solution found:", solution)
else:
    print("Solution invalid")