import math

def sum_prime_factors(n):
    """
    Calculates the sum of the unique prime factors of a given integer n.
    For example, for n=12 (2*2*3), the unique prime factors are 2 and 3, so the sum is 5.
    """
    if n < 2:
        return 0
    
    factors = set()
    d = 2
    temp_n = n
    
    # Handle factor of 2
    if temp_n % 2 == 0:
        factors.add(2)
        while temp_n % 2 == 0:
            temp_n //= 2
            
    # Handle odd factors
    d = 3
    while d * d <= temp_n:
        if temp_n % d == 0:
            factors.add(d)
            while temp_n % d == 0:
                temp_n //= d
        d += 2
        
    # If temp_n is still greater than 1, it must be a prime factor
    if temp_n > 1:
        factors.add(temp_n)
        
    return sum(factors)

def find_the_number():
    """
    Finds the smallest integer N that is the sum of two different integers a and b,
    where the sum of the prime factors of a and b is equal to 20.
    """
    # We are looking for the smallest N = a + b
    # Let's iterate through possible sums N, starting from the smallest possible sum (2+3=5)
    # The smallest integers a and b will likely be used, so N won't be extremely large.
    max_n_to_check = 100 
    
    for n_sum in range(4, max_n_to_check):
        # For each sum N, check pairs (a, b) such that a+b = N and a < b
        # Let's iterate 'a' from 2 up to n_sum / 2
        for a in range(2, (n_sum // 2) + 1):
            b = n_sum - a
            
            # Condition from not R: a and b must be different.
            # Our loop 'a' goes up to n_sum/2, so 'a' will be <= 'b'.
            # If we check a != b, it covers the condition.
            if a == b:
                continue

            spf_a = sum_prime_factors(a)
            spf_b = sum_prime_factors(b)

            # Condition from not Q: The sum of prime factors is 20.
            if spf_a + spf_b == 20:
                # Since we are iterating n_sum upwards, the first N we find is the smallest.
                # This satisfies condition not P.
                print(f"The number N that satisfies the negation is {n_sum}.")
                print(f"The two different integers are a = {a} and b = {b}.")
                print(f"The sum of prime factors for a={a} is {spf_a}.")
                print(f"The sum of prime factors for b={b} is {spf_b}.")
                print(f"The total sum of prime factors is {spf_a} + {spf_b} = 20.")
                print(f"The final equation for N is: {a} + {b} = {n_sum}")
                return n_sum

    return None

# Execute the function to find and print the solution
found_n = find_the_number()

# The final answer part is for the platform, the script prints the reasoning.
# if found_n is not None:
#    print(f"\n<<< {found_n} >>>")
