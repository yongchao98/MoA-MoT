def compute_euler_char_mod_k(k: int):
    """
    Computes the reduced Euler characteristic of the complex Delta_k modulo k.

    Args:
        k: A prime number such that k >= 3.
    """
    # First, we validate that k is a prime number >= 3.
    if not isinstance(k, int) or k < 3:
        print(f"Error: k={k} must be an integer greater than or equal to 3.")
        return
    is_prime = all(k % i != 0 for i in range(2, int(k**0.5) + 1))
    if not is_prime:
        print(f"Error: k={k} must be a prime number.")
        return

    # The reduced Euler characteristic of Delta_k modulo k is given by the formula (k-1)/2.
    # Since k is an odd prime, k-1 is even, so the division results in an integer.
    # Let R = (k-1)/2. Since 0 < R < k, R mod k is just R.
    
    # The equation for the result
    numerator = k - 1
    denominator = 2
    result = numerator // denominator

    print(f"For k = {k}:")
    print(f"The calculation for chi_hat(Delta_{k}) mod {k} is:")
    print(f"({k} - 1) / {denominator} = {numerator} / {denominator} = {result}")
    print(f"The final answer is {result}.")
    print("-" * 30)

if __name__ == '__main__':
    # We can test the function for several prime values of k.
    # The problem is posed for a general prime k, so we show a few examples.
    print("Computing the reduced Euler characteristic chi_hat(Delta_k) mod k for some primes k >= 3.")
    
    compute_euler_char_mod_k(3)
    compute_euler_char_mod_k(5)
    compute_euler_char_mod_k(7)
    compute_euler_char_mod_k(11)
    compute_euler_char_mod_k(13)
