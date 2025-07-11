def solve():
    """
    Calculates the next three capital letters in the sequence.
    """
    # The first letter of the next three terms is 'O'.
    # The numerical value for 'O' (A=0, B=1, ...) is 14.
    n1 = 14
    
    # The second letters are assumed to be 'A', 'B', 'C' for the next three terms.
    # Their numerical values are 0, 1, 2.
    next_n2_values = [0, 1, 2]
    
    results = []
    
    print("The formula is: L3 = (12 * L1 + 16 * L2 + 5) % 26")
    print("Where L1, L2, L3 are the 0-indexed values of the letters (A=0, Z=25).\n")
    
    for n2 in next_n2_values:
        # The derived formula
        n3 = (12 * n1 + 16 * n2 + 5) % 26
        
        # Convert numerical values back to letters
        l1 = chr(ord('A') + n1)
        l2 = chr(ord('A') + n2)
        l3 = chr(ord('A') + n3)
        
        # Print the calculation
        print(f"For the triplet starting with {l1}{l2}:")
        print(f"L1 = {n1}, L2 = {n2}")
        print(f"L3 = (12 * {n1} + 16 * {n2} + 5) % 26 = {n3}")
        print(f"Resulting triplet: {l1}{l2}{l3}\n")
        
        results.append(f"{l1}{l2}{l3}")
        
    print("The next three capital letters are:")
    print(" ".join(results))

solve()
<<<OAR OBH OCX>>>