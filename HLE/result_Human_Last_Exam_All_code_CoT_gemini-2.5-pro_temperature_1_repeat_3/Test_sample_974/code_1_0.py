import math

def get_prime_divisors(n):
    """
    Returns a list of unique prime divisors of a number n.
    """
    factors = set()
    # Check for 2
    if n % 2 == 0:
        factors.add(2)
        while n % 2 == 0:
            n //= 2
    # Check for odd factors
    d = 3
    while d * d <= n:
        if n % d == 0:
            factors.add(d)
            while n % d == 0:
                n //= d
        d += 2
    if n > 1:
        factors.add(n)
    return sorted(list(factors))

def solve_and_explain():
    """
    Solves the problem and prints the reasoning and the final answer.
    """
    q = 12740347

    # Step 1: Find all prime divisors of q.
    prime_divisors_of_q = get_prime_divisors(q)

    # Step 2: Explain the problem and the setup.
    print(f"Let q = {q}.")
    print("We need to find all prime divisors p of q such that the number of elements of order p is equal in PSL(3, q^2) and PSL(4, q).")
    
    # Step 3: State the relationship between N_p(PSL) and N_p(SL).
    print("\nThe number of elements of order p, N_p(G), in a group G is considered.")
    print("The condition is: N_p(PSL(3, q^2)) = N_p(PSL(4, q)).")
    print("The orders of the centers are |Z(SL(3, q^2))| = gcd(3, q^2 - 1) = 3 and |Z(SL(4, q))| = gcd(4, q - 1) = 2.")
    print("Any prime divisor p of q is not 2 or 3, so p is coprime to the orders of the centers.")
    print("A theorem states that if p does not divide |Z(G)|, then N_p(G/Z(G)) = N_p(G).")
    print("Thus, the condition simplifies to: N_p(SL(3, q^2)) = N_p(SL(4, q)).")
    
    # Step 4: Count the elements of order p.
    print("\nFor a prime p = char(k), elements of order p are non-identity unipotent matrices whose largest Jordan block size m satisfies 1 < m <= p.")
    print(f"Any prime divisor p of q={q} is much larger than 3 or 4.")
    print("Therefore, this condition holds for all non-identity unipotent matrices in both groups.")
    
    # Step 5: State the formulas for the number of elements.
    print("\nThe number of unipotent elements in SL(n, k) is k^(n*(n-1)).")
    print("For SL(3, q^2), the number of elements of order p is (q^2)^(3*(3-1)) - 1.")
    print("For SL(4, q), the number of elements of order p is q^(4*(4-1)) - 1.")

    # Step 6: Show the final equation.
    print("\nThe final equation is:")
    print(f"({q}^2)^(3*2) - 1 = {q}^(4*3) - 1")
    print(f"{q}^12 - 1 = {q}^12 - 1")
    print("This equality is always true for q > 1.")
    
    # Step 7: Conclude and print the answer.
    print("\nThis means the condition holds for ALL prime divisors of q.")
    print(f"The prime divisors of q = {q} are:")
    for p in prime_divisors_of_q:
        print(p)

solve_and_explain()
<<<12740347>>>