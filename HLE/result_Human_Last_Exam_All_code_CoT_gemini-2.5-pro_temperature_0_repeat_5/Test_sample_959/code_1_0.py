import math

def euler_totient(n):
    """Computes Euler's totient function phi(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    # For a prime number p, phi(p) = p - 1.
    # A simple primality test for illustration.
    is_prime = True
    if n <= 1:
        is_prime = False
    else:
        for i in range(2, int(math.sqrt(n)) + 1):
            if n % i == 0:
                is_prime = False
                break
    
    if is_prime:
        return n - 1
    
    # General case (not needed here but good for completeness)
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

# The group G in the problem is the trivial group {1}.
# The set of central extensions E of G by C=Z_31 has only one element, E isomorphic to Z_31.
# We need to compute the order of the outer automorphism group of E = Z_31.
# o(E) = |Out(Z_31)| = |Aut(Z_31)| since Z_31 is abelian.
# |Aut(Z_31)| = phi(31).

p = 31
order_of_outer_automorphism_group = euler_totient(p)

# The sum is over a single element E.
total_sum = order_of_outer_automorphism_group

print(f"The group G is trivial.")
print(f"There is only one central extension E, and E is isomorphic to the cyclic group C of order {p}.")
print(f"The order of the outer automorphism group of E is phi({p}).")
print(f"The final equation is: phi({p}) = {p} - 1 = {order_of_outer_automorphism_group}")
print(f"The sum is {total_sum}.")
