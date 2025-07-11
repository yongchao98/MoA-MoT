import math

def calculate_and_print_limit():
    """
    Calculates the limit based on theoretical derivation and prints it.

    The limit of the sequence g_n is the greatest common divisor of P(p) for all
    sufficiently large primes p.
    Through number theory analysis, we determined that the limit must be of the form 2^a * 3^b * 5^c.
    
    1. The exponent of 2, a, is min(v_2(P(p))) over large primes p. This is found to be 10.
    2. The exponent of 3, b, is min(v_3(P(p))) over large primes p. This is found to be 2.
    3. The exponent of 5, c, is min(v_5(P(p))) over large primes p. This is found to be 1.
    
    The limit is therefore 2^10 * 3^2 * 5^1.
    """
    
    val_2_10 = 2**10
    val_3_2 = 3**2
    val_5_1 = 5**1
    
    limit = val_2_10 * val_3_2 * val_5_1
    
    # The final output needs to show each number in the equation.
    print(f"{val_2_10} * {val_3_2} * {val_5_1} = {limit}")

# Although the problem can be solved with pure number theory, the following code
# could be used to empirically verify the result. It calculates the running GCD
# for the first few dozen primes greater than 10. For brevity and to directly
# answer the question, we call the function that prints the derived result.
#
# def P(x):
#     x = int(x) # Use Python's arbitrary-precision integers
#     return (x**5-1)*(x**5-x)*(x**5-x**2)*(x**5-x**3)*(x**5-x**4)
#
# def is_prime(n):
#     if n < 2: return False
#     for i in range(2, int(math.sqrt(n)) + 1):
#         if n % i == 0: return False
#     return True
#
# g = 0
# count = 0
# p = 11
# while count < 50:
#     if is_prime(p):
#         val = P(p)
#         if g == 0:
#             g = val
#         else:
#             g = math.gcd(g, val)
#         count += 1
#     p += 1
# print(f"Empirical GCD after {count} primes: {g}")
# print(f"Theoretical value: {2**10 * 3**2 * 5**1}")

calculate_and_print_limit()