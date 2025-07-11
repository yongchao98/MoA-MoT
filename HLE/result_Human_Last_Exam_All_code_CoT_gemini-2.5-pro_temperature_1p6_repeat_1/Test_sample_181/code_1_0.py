import math

def euler_totient(n):
    """Calculates Euler's totient function phi(n)."""
    # For a prime number p, phi(p) = p - 1.
    if n == 7:
        return n - 1
    # General implementation for completeness, though not strictly needed here.
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
    return int(result)

def solve_foo_7():
    """
    Calculates the value of foo(7) based on group theory principles.
    """
    n = 7

    # The function foo(n) counts the number of group structures on a labeled set of n elements.
    # For a prime n, there's only one group up to isomorphism: the cyclic group C_n.
    # The count is given by the formula: n! / |Aut(C_n)|.
    
    # |Aut(C_n)| is given by Euler's totient function, phi(n).
    aut_cn_size = euler_totient(n)
    
    # n! is n factorial.
    n_factorial = math.factorial(n)
    
    # Calculate the final result.
    result = n_factorial // aut_cn_size

    print("The function foo(7) computes the number of distinct group structures on a set of 7 elements.")
    print("This value is calculated using the formula: n! / |Aut(C_n)|, for n=7.")
    print(f"For n=7, n! is {n_factorial}.")
    print(f"The size of the automorphism group, |Aut(C_7)|, is phi(7) = {aut_cn_size}.")
    print("The final equation is:")
    print(f"{n_factorial} / {aut_cn_size} = {result}")

solve_foo_7()