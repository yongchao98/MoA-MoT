import math

def euler_totient(n):
    """
    Computes Euler's totient function, phi(n).
    This function counts the number of positive integers up to a given integer n
    that are relatively prime to n.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
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

# The problem concerns central extensions by C_31, the cyclic group of order 31.
n = 31

# As derived in the explanation, the group G is trivial.
# This implies there is only one central extension E, up to isomorphism, and E is isomorphic to C_31.
# The set of extensions E has only one element.
# We need to compute the order of the outer automorphism group of E, o(E).
# o(E) = |Out(E)| = |Aut(E)| since E = C_31 is abelian.
# The order of the automorphism group of C_n is given by Euler's totient function, phi(n).
# So, we need to calculate phi(31).
order_of_out_E = euler_totient(n)

# The sum is over a single element, so the sum is just this value.
final_sum = order_of_out_E

print("Step 1: The group G is the trivial group {1}.")
print("Step 2: There is a unique central extension E, and E is isomorphic to the cyclic group C_31.")
print("Step 3: The order of the outer automorphism group of E is phi(31).")
print(f"The calculation yields phi({n}) = {order_of_out_E}.")
print("\nThe final sum is over a single term.")
print(f"The equation for the sum is: Sum = o(E) = {final_sum}")
print(f"The computed sum is {final_sum}.")
