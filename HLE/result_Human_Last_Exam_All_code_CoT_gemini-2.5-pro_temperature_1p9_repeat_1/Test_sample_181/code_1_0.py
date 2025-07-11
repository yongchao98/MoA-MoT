import math

def solve():
    """
    This function calculates the value of foo(7).

    The provided function `foo(n)` is designed to count the number of valid group structures on a set of `n` elements.
    It iterates through all possible `n x n` multiplication tables and checks if they satisfy the group axioms:
    1. Associativity: (a*b)*c = a*(b*c)
    2. Identity Element: There exists a left identity 'e' such that e*x = x for all x.
    3. Inverse Element: For each element 'a', there exists a right inverse 'x' such that a*x = e.

    A structure satisfying these axioms is a group.

    Running foo(7) directly is computationally infeasible because it would have to check 7^(7*7) = 7^49 possible multiplication tables.

    Instead, we can solve this using group theory. The problem is equivalent to finding the number of distinct group tables of order 7.
    - Any group of a prime order 'p' is cyclic and isomorphic to Z_p. Since 7 is prime, any group of order 7 is isomorphic to Z_7.
    - The number of distinct group structures on n elements isomorphic to a group G is given by n! / |Aut(G)|, where |Aut(G)| is the size of the automorphism group of G.
    - For G = Z_n, the size of the automorphism group is given by Euler's totient function, phi(n).
    - For a prime number p, phi(p) = p - 1.
    - Thus, for n=7, |Aut(Z_7)| = phi(7) = 7 - 1 = 6.
    - The total number of group structures on 7 elements is 7! / phi(7).
    """
    n = 7

    # For a prime number n, phi(n) is n - 1.
    automorphisms_count = n - 1

    # Calculate n factorial.
    n_factorial = math.factorial(n)

    # The result is n! divided by the number of automorphisms.
    result = n_factorial // automorphisms_count

    print(f"{n_factorial} / {automorphisms_count} = {result}")

solve()