import math

def solve():
    """
    Solves the problem by determining the limit of the sequence g_n.
    The method is to find the p-adic valuation of the limit L for each prime p.
    """
    print("Let L be the limit of the sequence (g_n). We determine L by finding its prime factorization.")
    print("The exponent of a prime q in the factorization of L is v_q(L) = lim_{n->inf} min_{p>n} v_q(P(p)).\n")

    # Step 1: Analyze for primes q > 5
    print("Step 1: Analyzing prime factors q > 5.")
    print("A prime q can divide L only if P(a) is divisible by q for all integers a.")
    print("P(X) = (X^5-1)(X^5-X)(X^5-X^2)(X^5-X^3)(X^5-X^4)")
    print("P(a) is divisible by q if a^5 = a^k (mod q) for some k in {0,1,2,3,4}.")
    print("This is equivalent to a^(5-k) = 1 (mod q) for a != 0.")
    print("This means the order of 'a' in (Z/qZ)* must divide one of {1,2,3,4,5}.")
    print("For any prime q > 5, the group (Z/qZ)* is cyclic of order q-1 > 4.")
    print("Let 'a' be a generator. Its order is q-1, which does not divide any of {1,2,3,4,5}.")
    print("So, P(a) is not divisible by q.")
    print("By Dirichlet's theorem, there are infinitely many primes p = a (mod q).")
    print("For these primes, v_q(P(p)) = 0. Thus, min_{p>n} v_q(P(p)) = 0 for all n.")
    v_q_L_gt_5 = 0
    print(f"Therefore, v_q(L) = {v_q_L_gt_5} for all primes q > 5.\n")

    # Step 2: Analyze for q = 2
    print("Step 2: Analyzing v_2(L).")
    print("For any odd prime p, p^2 = 1 (mod 8), so p^5 = p (mod 8).")
    print("The valuations of the factors of P(p) are:")
    print("v_2(p^5-1) = v_2(p-1)")
    print("v_2(p^5-p) = v_2(p(p^4-1)) = v_2(p(p^2-1)(p^2+1)) = v_2(p^2-1) + v_2(p^2+1) = 3 + 1 = 4")
    print("v_2(p^5-p^2) = v_2(p^2(p^3-1)) = v_2(p-1)")
    print("v_2(p^5-p^3) = v_2(p^3(p^2-1)) = 3")
    print("v_2(p^5-p^4) = v_2(p^4(p-1)) = v_2(p-1)")
    print("v_2(P(p)) = v_2(p-1) + 4 + v_2(p-1) + 3 + v_2(p-1) = 7 + 3*v_2(p-1).")
    print("The minimum value occurs when v_2(p-1) is minimal, which is 1 (for p = 3 mod 4).")
    v_2_L = 7 + 3 * 1
    print(f"The minimum value is 7 + 3*1 = {v_2_L}.")
    print(f"By Dirichlet's theorem, such primes exist for any large n. So, v_2(L) = {v_2_L}.\n")

    # Step 3: Analyze for q = 3
    print("Step 3: Analyzing v_3(L).")
    print("For a prime p != 3:")
    print("Case p = 1 (mod 3): v_3(P(p)) = 5*v_3(p-1) + 1. Minimum value (for v_3(p-1)=1) is 6.")
    print("Case p = 2 (mod 3): v_3(P(p)) = 2*v_3(p+1). Minimum value (for v_3(p+1)=1) is 2.")
    v_3_L = 2
    print(f"The overall minimum is {v_3_L}, achieved for primes p = 2, 5 (mod 9).")
    print(f"By Dirichlet's theorem, such primes exist for any large n. So, v_3(L) = {v_3_L}.\n")

    # Step 4: Analyze for q = 5
    print("Step 4: Analyzing v_5(L).")
    print("For a prime p != 5:")
    print("Case p = 1 (mod 5): v_5(P(p)) = 5*v_5(p-1) + 1. Minimum value (for v_5(p-1)=1) is 6.")
    print("Case p != 1 (mod 5): Only the factor (p^5-p) is divisible by 5. v_5(P(p)) = v_5(p^4-1).")
    print("The minimum of v_5(p^4-1) is 1 (e.g., for p = 2 mod 25).")
    v_5_L = 1
    print(f"The overall minimum is {v_5_L}.")
    print(f"By Dirichlet's theorem, such primes exist for any large n. So, v_5(L) = {v_5_L}.\n")

    # Step 5: Calculate the final limit L
    print("Step 5: Calculating the final limit L.")
    p2 = 2**v_2_L
    p3 = 3**v_3_L
    p5 = 5**v_5_L
    limit = p2 * p3 * p5
    
    print(f"The limit L is the product of these prime powers:")
    print(f"L = 2^{v_2_L} * 3^{v_3_L} * 5^{v_5_L}")
    print(f"L = {p2} * {p3} * {p5}")
    print(f"L = {limit}")

solve()