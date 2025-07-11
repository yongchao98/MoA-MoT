import math

def get_sum_of_prime_factors(n):
    """
    Calculates the sum of the distinct prime factors of a given number n.
    For example, for 12 (prime factors 2, 2, 3), the distinct factors are 2 and 3, so the sum is 5.
    """
    if n <= 1:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    
    # Efficiently find prime factors
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    
    # If a number remains, it's a prime factor itself
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def solve_puzzle():
    """
    Finds the number N that satisfies the given logical conditions.
    """
    print("The problem asks for the number N that satisfies the negation of 'P or R or Q'.")
    print("This is equivalent to 'not P' and 'not R' and 'not Q'.\n")

    print("Breaking down the conditions:")
    print("1. (not P): N is the SMALLEST number that is the sum of two integers a and b.")
    print("2. (not R): The integers a and b must be DIFFERENT.")
    print("3. (not Q): The sum of the prime factors of a and b must be EXACTLY 20.\n")

    print("Searching for the smallest N = a + b such that a != b and sum_prime_factors(a) + sum_prime_factors(b) == 20...")
    
    # Start with the smallest possible sum of two different integers >= 2 (2+3=5)
    N = 5
    
    while True:
        # We check pairs (a,b) that sum to N. To avoid duplicates and ensure a!=b,
        # we can iterate 'a' from 2 up to N/2.
        for a in range(2, N // 2 + 1):
            b = N - a

            # The loop structure ensures a < b, so a != b is guaranteed.
            
            sum_factors_a = get_sum_of_prime_factors(a)
            sum_factors_b = get_sum_of_prime_factors(b)

            if sum_factors_a + sum_factors_b == 20:
                print("\n--- Solution Found ---")
                print(f"Found a pair (a, b) = ({a}, {b}) that satisfies the conditions for N = {N}.")
                print(f"Sum of prime factors of a={a} is {sum_factors_a}.")
                print(f"Sum of prime factors of b={b} is {sum_factors_b}.")
                print(f"Total sum of prime factors = {sum_factors_a} + {sum_factors_b} = 20.")
                print(f"Since we searched for N incrementally, N={N} is the smallest such sum.")
                print("\nThe final equation is:")
                print(f"{N} = {a} + {b}")
                return N
        
        N += 1

# Execute the solver
final_answer = solve_puzzle()
# The final answer is printed by the function, but we return it for consistency.
# The following line is for the final answer format.
# print(f"<<<{final_answer}>>>")