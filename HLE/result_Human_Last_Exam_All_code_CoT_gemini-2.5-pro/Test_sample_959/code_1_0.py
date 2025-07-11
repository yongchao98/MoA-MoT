import math

def euler_totient(n):
    """Computes Euler's totient function phi(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer.")
    
    # Since n=31 is a prime number, phi(n) = n - 1.
    # For a general implementation:
    count = 0
    for i in range(1, n + 1):
        if math.gcd(n, i) == 1:
            count += 1
    return count

# The problem simplifies to calculating the order of the automorphism group of C_31.
# G is the trivial group {1}.
# The set of central extensions E of G by C=C_31 has only one element, E isomorphic to C_31.
# The order of the outer automorphism group of E, o(E), is |Aut(C_31)| since C_31 is abelian.
# |Aut(C_31)| is given by Euler's totient function phi(31).
p = 31
order = euler_totient(p)

# The sum is over a single element, so the sum is just the order itself.
# The final equation is: Sum = phi(31) = 30
print("The group G is the trivial group {1}.")
print("There is only one central extension E, which is isomorphic to C_31.")
print("The order of the outer automorphism group of E is phi(31).")
print(f"phi({p}) = {p-1}")
print(f"The sum is therefore {order}.")
