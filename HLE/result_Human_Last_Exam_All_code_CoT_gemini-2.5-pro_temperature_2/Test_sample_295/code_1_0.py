def is_prime(n):
    """A helper function to check if a number is prime."""
    if not isinstance(n, int) or n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    # Check for odd divisors from 3 up to sqrt(n)
    for i in range(3, int(n**0.5) + 1, 2):
        if n % i == 0:
            return False
    return True

def compute_chi_mod_k(k):
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k modulo k.

    Args:
        k (int): A prime number such that k >= 3.
    """
    if not (isinstance(k, int) and k >= 3 and is_prime(k)):
        print(f"Error: The input k must be a prime number greater than or equal to 3. Received: {k}")
        return

    # Based on the mathematical derivation, the result is given by the formula (k-1)/2.
    k_minus_1 = k - 1
    result = k_minus_1 // 2

    # As requested, we print the numbers in the final equation.
    print(f"For the prime k = {k}:")
    print("The reduced Euler characteristic of Delta_k modulo k is calculated with the formula:")
    print(f"  \u03C7\u0302(\u0394_{k}) mod k = (k - 1) / 2")
    # Displaying the actual calculation
    print("Substituting the value of k:")
    print(f"  \u03C7\u0302(\u0394_{{{k}}}) mod {k} = ({k} - 1) / 2 = {k_minus_1} / 2 = {result}")
    
    # The result of `result mod k` is just `result` itself, as 0 < result < k.
    print(f"The final result is {result}.")
    print("-" * 30)

if __name__ == '__main__':
    # The problem is posed for a general prime k >= 3.
    # We demonstrate the solution for a few such prime numbers.
    example_primes = [3, 5, 7, 11, 17]

    for prime_k in example_primes:
        compute_chi_mod_k(prime_k)
