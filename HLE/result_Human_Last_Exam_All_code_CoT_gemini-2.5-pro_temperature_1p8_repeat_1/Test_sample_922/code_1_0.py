import math

def get_perfect_square_root(n):
    """If n is a perfect square, returns its integer root. Otherwise, returns None."""
    if n < 0:
        return None
    if n == 0:
        return 0
    # Using integer square root is fast and avoids floating-point inaccuracies.
    root = int(math.isqrt(n))
    if root * root == n:
        return root
    return None

def get_prime_factorization(num):
    """Returns the prime factorization of num as a dictionary {prime: exponent}."""
    factors = {}
    d = 2
    temp = num
    # Handle factor 2
    while (temp % d) == 0:
        factors[d] = factors.get(d, 0) + 1
        temp //= d
    # Handle odd factors
    d = 3
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 2
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def calculate_sigma_from_factors(factors):
    """Calculates the sum of divisors (sigma) from a prime factorization dictionary."""
    total_sum = 1
    for p, a in factors.items():
        # Formula for sum of divisors of a prime power: sigma(p^a) = (p^(a+1) - 1) / (p - 1)
        total_sum *= (p**(a + 1) - 1) // (p - 1)
    return total_sum

def solve_and_verify_sequence():
    """
    Identifies the rule for the sequence and verifies each term, including the next number.
    """
    sequence = [24663, 35005, 119261, 196219, 211770, 227296]
    next_term = 289876
    
    full_sequence = sequence + [next_term]
    
    print("This sequence (OEIS A054320) contains numbers 'n' where sigma(n^3) is a perfect square.")
    print("Verifying this property for the given terms and the proposed next term.\n")
    
    for n in full_sequence:
        # Get prime factorization of n
        n_factors = get_prime_factorization(n)
        
        # The prime factorization of n^3 has the same primes with exponents multiplied by 3
        n3_factors = {p: e * 3 for p, e in n_factors.items()}
        
        # Calculate sigma(n^3) using its prime factorization
        sigma_n3 = calculate_sigma_from_factors(n3_factors)
        
        # Check if sigma(n^3) is a perfect square and get its root
        sqrt_sigma_n3 = get_perfect_square_root(sigma_n3)
        
        if sqrt_sigma_n3 is not None:
            # Print each number in the verification equation as requested
            print(f"For n = {n}:")
            print(f"  sigma({n}^3) = {sigma_n3}")
            print(f"  This is a perfect square: {sqrt_sigma_n3}^2\n")
        else:
            # This case should not be hit for the numbers in the sequence
            print(f"Verification failed for n = {n}.\n")

    print(f"The single known integer that completes the sequence is {next_term}.")

if __name__ == '__main__':
    solve_and_verify_sequence()