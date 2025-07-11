import math

def phi(n):
    """
    Computes Euler's totient function. For a prime p, phi(p) = p - 1.
    """
    if n < 1:
        return 0
    
    # Check if n is prime for a simple calculation
    is_prime = True
    if n == 1:
        is_prime = False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            is_prime = False
            break
            
    if is_prime:
        return n - 1
    
    # General case (not strictly needed for n=31)
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

# The order of the cyclic group C
n = 31

# As explained in the steps above, the group G is trivial.
# This simplifies the problem to finding the order of the outer automorphism group
# of the single central extension E = C = Z_n.

# o(E) = |Out(E)| = |Aut(E)| / |Inn(E)|
# For E = Z_n:
# |Aut(Z_n)| is given by Euler's totient function, phi(n).
order_of_Aut_C = phi(n)

# |Inn(Z_n)| is 1 because Z_n is abelian.
order_of_Inn_C = 1

# Calculate the order of the outer automorphism group.
order_of_Out_C = order_of_Aut_C // order_of_Inn_C

# The collection E has only one element, so the sum is just this value.
total_sum = order_of_Out_C

print("The problem asks for the sum of the orders of outer automorphism groups for all central extensions E.")
print("Based on the analysis, the group G is trivial, and there is only one such extension E, which is isomorphic to C = Z_31.")
print("The sum is therefore just the order of Out(Z_31).")
print("\nThe final equation is:")
print(f"Sum = |Out(Z_{n})|")
print(f"|Out(Z_{n})| = |Aut(Z_{n})| / |Inn(Z_{n})|")
print(f"|Aut(Z_{n})| = phi({n}) = {order_of_Aut_C}")
print(f"|Inn(Z_{n})| = {order_of_Inn_C}")
print(f"|Out(Z_{n})| = {order_of_Aut_C} / {order_of_Inn_C} = {order_of_Out_C}")
print(f"Sum = {total_sum}")