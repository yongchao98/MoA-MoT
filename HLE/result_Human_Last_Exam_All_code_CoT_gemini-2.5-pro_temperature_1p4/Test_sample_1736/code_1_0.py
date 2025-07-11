import math

def get_sum_of_unique_prime_factors(n):
    """
    Calculates the sum of unique prime factors of a number n.
    For example, for 12 (2*2*3), the unique prime factors are {2, 3}, so the sum is 2+3=5.
    """
    if n <= 1:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 1
    if temp_n > 1:
        factors.add(temp_n)
    return sum(factors)

def find_the_number():
    """
    Finds the number N based on the logical conditions derived from the problem.
    ¬(P ∨ R ∨ Q) <=> ¬P ∧ ¬R ∧ ¬Q
    - ¬R: a != b
    - ¬Q: sum_prime_factors(a) + sum_prime_factors(b) = 20
    - ¬P: N = a + b is the smallest such sum.
    """
    target_sum = 20
    
    # We iterate through possible sums N, starting from the smallest (1+2=3)
    N = 3
    while True:
        # For each sum N, check all pairs (a, b) where a+b=N and a < b
        for a in range(1, math.ceil(N / 2)):
            b = N - a
            
            # This check ensures a and b are different
            if a == b:
                continue

            sum_factors_a = get_sum_of_unique_prime_factors(a)
            sum_factors_b = get_sum_of_unique_prime_factors(b)
            
            # This is condition ¬Q
            if sum_factors_a + sum_factors_b == target_sum:
                # Because we are iterating N upwards, the first N we find is the smallest.
                # This satisfies condition ¬P.
                print(f"Found the smallest number N that satisfies the conditions.")
                print(f"The number is the sum of a and b, where a={a} and b={b}.")
                print("\nVerification:")
                print(f"- Condition ¬R (a != b): {a} != {b} is True.")
                print(f"- Condition ¬Q (sum of prime factors is 20):")
                print(f"  - Sum of prime factors of {a}: {sum_factors_a}")
                print(f"  - Sum of prime factors of {b}: {sum_factors_b}")
                print(f"  - Total sum: {sum_factors_a} + {sum_factors_b} = {sum_factors_a + sum_factors_b}")
                print(f"- Condition ¬P (N is the smallest sum): N = {N} is the first value found by ascending search.")
                
                print(f"\nThe final equation is {N} = {a} + {b}")
                return N
                
        N += 1

# Execute the search and print the final result
final_N = find_the_number()
print(f"<<<{final_N}>>>")