import math

def get_divisors(n):
    """Returns a list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return list(divs)

def sigma(n, k):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if not isinstance(n, int) or n <= 0:
        raise ValueError("Input n must be a positive integer.")
    divs = get_divisors(n)
    return sum(d**k for d in divs)

def get_f_coeff(n):
    """Computes the n-th coefficient of the cusp form f."""
    if n % 2 == 1:
        # n is odd
        return sigma(n, 7)
    else:
        # n is even
        return sigma(n, 7) - sigma(n // 2, 7)

def main():
    """
    Calculates the sum of the first three non-zero coefficients
    of the specified cusp form f.
    """
    # The coefficients are guaranteed to be non-zero for n > 0.
    # a_n = sigma_7(n) for n odd
    # a_n = sigma_7(n) - sigma_7(n/2) for n even

    # Calculate the first three coefficients
    a1 = get_f_coeff(1)
    a2 = get_f_coeff(2)
    a3 = get_f_coeff(3)

    # Calculate the sum
    total_sum = a1 + a2 + a3

    # Print the result in the requested format
    print(f"The first three non-zero coefficients are a_1, a_2, and a_3.")
    print(f"a_1 = sigma_7(1) = {a1}")
    print(f"a_2 = sigma_7(2) - sigma_7(1) = (1^7 + 2^7) - 1^7 = 129 - 1 = {a2}")
    print(f"a_3 = sigma_7(3) = 1^7 + 3^7 = 1 + 2187 = {a3}")
    print(f"The sum is: {a1} + {a2} + {a3} = {total_sum}")

if __name__ == "__main__":
    main()
