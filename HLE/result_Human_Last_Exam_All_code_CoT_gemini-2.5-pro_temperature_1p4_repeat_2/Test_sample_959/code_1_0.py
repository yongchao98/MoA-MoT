import math

def phi(n):
    """
    Computes Euler's totient function.
    For a prime number p, phi(p) = p - 1.
    """
    if n <= 0:
        return 0
    # A simple implementation for demonstration.
    # For a prime number, the calculation is trivial.
    # We can add a primality test for robustness, but for n=31 it's not necessary.
    is_prime = all(n % i != 0 for i in range(2, int(math.sqrt(n)) + 1))
    if is_prime:
        return n - 1
    
    # General case (not needed here but good practice)
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

def solve():
    """
    Solves the problem based on the reasoning that G is the trivial group.
    """
    # The group C is the cyclic group of order 31.
    n = 31

    # As explained, G is the trivial group.
    # This means there is only one central extension E up to isomorphism, and E is isomorphic to C.
    # We need to compute the order of the outer automorphism group of E = Z_31.
    # o(E) = |Out(Z_31)| = |Aut(Z_31)| since Z_31 is abelian.
    # |Aut(Z_31)| is given by Euler's totient function phi(31).
    
    # The final equation is: Sum = o(E) = phi(31)
    
    order_of_c = n
    order_of_outer_automorphism_group = phi(order_of_c)

    # The sum is over a single element E, so the sum is just o(E).
    total_sum = order_of_outer_automorphism_group
    
    print(f"The group G is the trivial group.")
    print(f"There is only one central extension E, which is isomorphic to the cyclic group C of order {order_of_c}.")
    print(f"The order of the outer automorphism group of E is phi({order_of_c}).")
    print(f"phi({order_of_c}) = {order_of_c} - 1 = {order_of_outer_automorphism_group}")
    print(f"The sum is therefore {total_sum}.")

solve()
