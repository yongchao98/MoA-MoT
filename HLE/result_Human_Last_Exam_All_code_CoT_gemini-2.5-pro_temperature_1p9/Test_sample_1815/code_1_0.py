import sympy

def find_primes_in_progression(a, n, limit):
    """Finds primes p <= limit such that p = a (mod n)."""
    primes = []
    # sympy.primerange is efficient for finding primes
    for p in sympy.primerange(1, limit):
        if (p % n) == a:
            primes.append(p)
    return primes

def check_property_q_for_set_A(n, limit):
    """
    Checks if the intersection of A = {p-1 | p is prime} and nZ is infinite.
    We demonstrate this by counting how many we find up to a certain limit.
    According to Dirichlet's Theorem, this number should grow without bound as the limit increases.
    """
    # A_intersect_nZ corresponds to primes p where p-1 is a multiple of n,
    # which is equivalent to p = 1 (mod n).
    primes_p_equiv_1_mod_n = find_primes_in_progression(1, n, limit)
    
    # The set A intersect nZ will have the same size as the set of these primes.
    intersection_size = len(primes_p_equiv_1_mod_n)
    
    print(f"For n = {n}:")
    print(f"Searching for primes p <= {limit} such that p = 1 (mod {n})...")
    print(f"Found {intersection_size} such primes.")
    if intersection_size > 10:
        print(f"Examples: {primes_p_equiv_1_mod_n[:10]}...")
    else:
        print(f"Examples: {primes_p_equiv_1_mod_n}")
    print("Dirichlet's Theorem states that this set is infinite.")
    print("-" * 20)

def main():
    """
    This program demonstrates that Property Q, a necessary condition for the existence
    of the described topologies, is false.
    Property Q: For any infinite set A of non-zero integers, there exists an n>0
    such that A contains only a finite number of multiples of n.

    We test this with the infinite set A = {p - 1 | p is a prime}.
    An element p-1 is a multiple of n if and only if p is congruent to 1 modulo n.
    Dirichlet's Theorem on Arithmetic Progressions states that there are infinitely
    many such primes for any n > 1.

    Therefore, the set (A intersect nZ) is infinite for all n > 1.
    This contradicts Property Q. Because Property Q is false, no such topology exists.
    The number of such topologies is 0.
    """
    print(main.__doc__)
    print("\nDemonstration for A = {p - 1 | p is prime}:")
    
    test_ns = [2, 3, 5, 10]
    limit = 1000
    for n in test_ns:
        check_property_q_for_set_A(n, limit)
        
    print("Since for any n > 1 we can find infinitely many elements in A that are multiples of n,")
    print("Property Q fails. This contradiction implies the number of such topologies is 0.")
    print("\nThe final answer is 0.")

if __name__ == '__main__':
    main()
