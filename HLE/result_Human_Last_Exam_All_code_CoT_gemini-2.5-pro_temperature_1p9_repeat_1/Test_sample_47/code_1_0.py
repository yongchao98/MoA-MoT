def solve_expected_time():
    """
    Calculates the expected time until a sequence appears from a random source.

    The method is based on identifying overlaps between prefixes and suffixes
    of the target sequence. The expected time E is the sum of n^k for all k
    where the prefix of length k matches the suffix of length k.
    """
    # Define the alphabet size and the target sequence
    n = 26
    S = "TENETENET"
    L = len(S)

    # Variables to store the calculation
    total_expected_time = 0
    equation_parts = []

    # Find all k where the prefix of length k matches the suffix of length k
    for k in range(1, L + 1):
        # Slicing in Python: S[:k] is the prefix of length k
        # S[-k:] is the suffix of length k
        if S[:k] == S[-k:]:
            term = f"{n}^{k}"
            equation_parts.append(term)
            total_expected_time += (n**k)
    
    # Format the final output
    equation = "E = " + " + ".join(equation_parts)

    print(f"The target sequence is: \"{S}\"")
    print(f"The alphabet size is: {n}")
    print("\nThe expected time (E) is calculated by finding all overlaps")
    print("where a prefix of the sequence is also a suffix.\n")
    print("The final equation is:")
    print(equation)
    print("\nCalculating the total expected time:")
    print(f"E = {total_expected_time:,}")

# Run the solver
solve_expected_time()

# The final result
final_answer = 26**9 + 26**7 + 26**5 + 26**1
# The following line is for the final answer extraction and is not printed to the user.
# print(f"<<<{final_answer}>>>")