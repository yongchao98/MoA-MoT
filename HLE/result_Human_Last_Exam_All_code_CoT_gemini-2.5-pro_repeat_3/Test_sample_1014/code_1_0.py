import math

def find_prime_factorization(num):
    """
    Finds the prime factors of a given number, prints the factorization equation,
    and returns the largest prime factor.
    """
    original_num = num
    factors = []

    # Handle the factor of 2
    while num % 2 == 0:
        factors.append(2)
        num = num // 2

    # Handle odd factors starting from 3
    # We only need to check up to the square root of the number
    for i in range(3, int(math.sqrt(num)) + 1, 2):
        while num % i == 0:
            factors.append(i)
            num = num // i

    # This condition is to handle the case when the remaining number
    # is a prime number greater than 2 (e.g., the last and largest factor)
    if num > 2:
        factors.append(num)

    # Print the factorization equation
    if factors:
        equation = f"{original_num} = {' * '.join(map(str, factors))}"
        print(equation)
        largest_prime = factors[-1]
        print(f"The largest prime divisor is: {largest_prime}")
    else:
        # This case handles if the original number itself is prime
        print(f"{original_num} is a prime number.")
        largest_prime = original_num
        print(f"The largest prime divisor is: {largest_prime}")
    
    return largest_prime

# The number to be factored
number_to_factor = 8139881
find_prime_factorization(number_to_factor)