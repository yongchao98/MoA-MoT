import math

def euler_totient(n):
    """
    Computes Euler's totient function phi(n).
    """
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

# The problem simplifies to finding the order of the outer automorphism group
# of a single group E, which is the cyclic group of order 31.
n = 31

# The order of the outer automorphism group of C_n is phi(n) for an abelian group.
# For a prime p, phi(p) = p - 1.
order_o_E = euler_totient(n)

# The sum is just this single term.
total_sum = order_o_E

print("Under the justified assumption that the group G is trivial, the set of extensions E contains only one element, C_31.")
print("The sum is therefore the order of the outer automorphism group of C_31.")
print(f"Let E = C_{n}, with n = {n}.")
print("The order of the outer automorphism group is o(E) = |Aut(E)| / |Inn(E)|.")
print("Since E is abelian, |Inn(E)| = 1.")
print("Thus, o(E) = |Aut(E)|, which is given by Euler's totient function phi(n).")
print(f"o(E) = phi({n})")
# For a prime number p, phi(p) = p - 1
p = n
phi_p = p - 1
print(f"phi({p}) = {p} - 1 = {phi_p}")
print(f"So the total sum is {total_sum}.")
