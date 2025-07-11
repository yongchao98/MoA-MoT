import math

def sum_of_prime_factors(n):
    """
    Calculates the sum of distinct prime factors of an integer n.
    For example, for 12 = 2*2*3, the distinct prime factors are 2 and 3, and the sum is 5.
    """
    if n < 2:
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

def solve():
    """
    Finds the smallest integer N that satisfies the given logical conditions.
    N must be the sum of two different integers a and b, where the sum 
    of the prime factors of a and b is exactly 20.
    """
    N = 4
    while True:
        # Iterate through pairs (a,b) where a + b = N and a < b.
        # We start 'a' from 2, as numbers less than 2 have no prime factors.
        for a in range(2, N // 2 + 1):
            if a == N - a:
                continue # a and b must be different.
            
            b = N - a
            
            sum_factors_a = sum_of_prime_factors(a)
            sum_factors_b = sum_of_prime_factors(b)
            
            # Check if the sum of prime factors is 20
            if sum_factors_a + sum_factors_b == 20:
                # Since we are iterating N upwards, the first solution found is the smallest.
                print(f"The first integer a is: {a}")
                print(f"The second integer b is: {b}")
                print(f"The smallest sum N is: {N}")
                # As per the instruction, printing the numbers in the final equation.
                print(f"{a} + {b} = {N}")
                return

        N += 1

if __name__ == '__main__':
    solve()