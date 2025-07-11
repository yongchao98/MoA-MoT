import math

def euler_totient(n):
    """
    Computes Euler's totient function phi(n).
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer")
    
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    return result

# The problem states that C is the cyclic group of order 31.
n = 31

# As argued in the text, there is only one central extension E in E,
# and E is isomorphic to C_31.
# The order of the outer automorphism group of E is phi(31).
o_E = euler_totient(n)

# The sum is over a single element, so it's just o(E).
total_sum = o_E

# As requested, printing the equation. Since there is only one element in E,
# the sum is just o(E).
print(f"o(E_1) = {o_E}")
print(f"Sum = {total_sum}")
