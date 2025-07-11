def phi(n):
    """
    Computes Euler's totient function phi(n), which counts the positive integers
    up to a given integer n that are relatively prime to n.
    For a prime number p, phi(p) = p - 1.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer")
    
    # A simple check for primality, sufficient for this problem's scale.
    is_prime = n > 1 and all(n % i != 0 for i in range(2, int(n**0.5) + 1))
    
    if is_prime:
        return n - 1
        
    # A general implementation for non-prime numbers.
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

def solve_group_problem():
    """
    Solves the user's problem by walking through the logical steps and calculating the final answer.
    """
    # The group C is the cyclic group of order 31.
    p = 31

    print("The solution follows from a mathematical analysis of the group G.")
    print("The main steps of the argument are:")
    print("1. The group G is shown to be the trivial group {1}.")
    print("   This is proven by showing that the relations force the order of each generator to be 1.")
    print("2. The set E of central extensions of G by C is classified by the second cohomology group H^2(G, C).")
    print(f"   Since G is trivial, H^2({{1}}, C) is also trivial.")
    print("3. This implies there is only one central extension E in the collection E, up to isomorphism.")
    print(f"   This unique extension E is isomorphic to C itself, so E is Z_{p}.")
    print("4. The sum collapses to a single term: the order of the outer automorphism group of E.")
    
    print("\nCalculating the final sum:")
    print(f"Sum = o(E) = |Out(Z_{p})|")

    # For an abelian group, the inner automorphism group is trivial.
    # So, Out(Z_p) is isomorphic to Aut(Z_p).
    print(f"|Out(Z_{p})| = |Aut(Z_{p})|")
    
    # The order of the automorphism group of a cyclic group Z_n is given by phi(n).
    print(f"|Aut(Z_{p})| = phi({p})")
    
    # Calculate phi(p)
    order_aut_group = phi(p)

    print(f"The equation for the final answer is: phi({p}) = {p} - 1")
    print(f"The result of the calculation is: {order_aut_group}")

    print("\n--------------------------------")
    print(f"The final sum is {order_aut_group}.")
    print("--------------------------------")


solve_group_problem()
<<<30>>>