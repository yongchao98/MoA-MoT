import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return result

# According to the reasoning, the problem reduces to calculating the order
# of the automorphism group of the cyclic group C of order 31.
# This is given by Euler's totient function, phi(31).
order_of_C = 31

# Calculate the order of the outer automorphism group of E.
# Since E is isomorphic to C_31 (an abelian group), Out(E) is isomorphic
# to Aut(E), and |Aut(C_31)| = phi(31).
order_of_outer_automorphism_group = phi(order_of_C)

# The set of extensions E has only one element, so the sum is just this value.
total_sum = order_of_outer_automorphism_group

print(f"The order of the cyclic group C is {order_of_C}.")
print(f"The order of the outer automorphism group o(E) is phi({order_of_C}).")
print(f"phi({order_of_C}) = {order_of_C} - 1 = {order_of_outer_automorphism_group}")
print(f"The final sum is {total_sum}.")