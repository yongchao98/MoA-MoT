def check_numbers(n1, n2):
    # Check all conditions for numbers
    if n1 in [8, 0, 9] or n2 in [8, 0, 9]:  # Guess 1 and 5
        return False
    if n1 < 4 or n2 < 4:  # Guess 3
        return False
    
    # Check guess 2 (50)
    if not ((n1 == 5 and n2 != 0) or (n2 == 5 and n1 != 0)):
        return False
    
    # Check guess 4 (76)
    if not ((n1 == 7 and n2 != 6) or (n2 == 7 and n1 != 6)):
        return False
    
    return True

def find_solution():
    # We know Z and V are the letters and they're in correct positions
    letters = ['Z', 'V']
    
    # Try all possible number combinations
    for n1 in range(10):
        for n2 in range(10):
            if n1 != n2 and check_numbers(n1, n2):
                print(f"Possible solution: [{n1}, {n2}, 'Z', 'V']")

find_solution()