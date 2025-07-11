import math

def find_largest_prime_divisor():
    """
    This function finds the largest prime divisor of a given number,
    prints the prime factorization, and returns the largest prime factor.
    """
    n = 8139881
    original_n = n
    factors = []

    # Handle the factor 2. For an odd number, this loop is skipped.
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Handle odd factors. Iterate from 3 up to the square root of n.
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n = n // i
        i += 2

    # This condition handles the case where n is a prime number
    # greater than 2 after all other divisions.
    if n > 2:
        factors.append(n)

    # Determine the largest prime divisor
    if not factors:
        # This case handles input n=1
        print("1 has no prime factors.")
        return None
    else:
        largest_prime = factors[-1]

    # Print the final equation and the result
    factor_string = " * ".join(map(str, factors))
    print(f"The prime factorization of {original_n} is: {original_n} = {factor_string}")
    print(f"The largest prime divisor of {original_n} is: {largest_prime}")
    
    return largest_prime

# Run the function and store the result for the final answer
largest_divisor = find_largest_prime_divisor()
# The final answer is formatted as requested.
# To avoid printing it here, we will place it after the code block.

# <<<8139881>>>