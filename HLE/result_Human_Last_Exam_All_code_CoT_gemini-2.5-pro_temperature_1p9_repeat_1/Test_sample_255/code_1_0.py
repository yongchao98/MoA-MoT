import math

def totient(n):
    """Computes Euler's totient function phi(n)."""
    if not isinstance(n, int) or n < 1:
        raise ValueError("Input must be a positive integer")
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

def get_divisors(n):
    """Computes divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve():
    """
    Solves for the dimension of the cohomology group H^2(G,M).
    The solution follows the Lyndon-Hochschild-Serre spectral sequence for the central extension
    1 -> Z -> G -> Q -> 1, where Z is the infinite cyclic group and Q = C_8 * C_8.
    """
    
    # Step 1: Compute dimensions of H^q(Z,M).
    # The action of Z = <z = a^8> on M is via T^8, where T is the permutation matrix.
    # M decomposes as a sum of Q[x]/Phi_d(x) for d|128.
    # H^q(Z, M) depends on T^8 - I. This operator is 0 on components M_d where d|8,
    # and an isomorphism on others.
    # So, H^0(Z, M) and H^1(Z, M) are isomorphic to the direct sum of M_d for d|8.
    
    divisors_of_8 = get_divisors(8)
    dim_V = sum(totient(d) for d in divisors_of_8)
    
    print("Step 1: Compute dimensions of H^q(Z,M). Let V be the resulting vector space.")
    print("Both H^0(Z,M) and H^1(Z,M) are isomorphic to a vector space V.")
    dim_eq = " + ".join([f"phi({d})" for d in divisors_of_8])
    dim_vals = " + ".join([f"{totient(d)}" for d in divisors_of_8])
    print(f"dim(V) = {dim_eq} = {dim_vals} = {dim_V}")
    print("dim H^q(Z,M) = 0 for q >= 2.")
    print("-" * 30)
    
    # The generators of Q act on V via T.
    # We compute H^p(C_8, V) to use in the Mayer-Vietoris sequence for Q = C_8 * C_8.
    
    # Step 2: Compute dimensions of H^p(C_8, V).
    # dim H^0(C_8, V) = dim ker(T-I) on V. T-I is zero only on the M_1 component.
    dim_H0_C8_V = totient(1)
    
    # dim H^{odd}(C_8,V) = dim(ker N / im(T-I)) = 0
    dim_ker_N = totient(2) + totient(4) + totient(8)
    dim_im_T_minus_I = totient(2) + totient(4) + totient(8)
    dim_H_odd_C8_V = dim_ker_N - dim_im_T_minus_I

    # dim H^{even>0}(C_8,V) = dim(ker(T-I) / im N) = 0
    dim_ker_T_minus_I = totient(1)
    dim_im_N = totient(1)
    dim_H_even_C8_V = dim_ker_T_minus_I - dim_im_N

    print("Step 2: Compute dimensions of H^p(C_8, V).")
    print(f"dim H^0(C_8, V) = phi(1) = {dim_H0_C8_V}")
    print(f"dim H^p(C_8, V) = 0 for p >= 1.")
    print("-" * 30)
    
    # Step 3: Compute dimensions of H^p(Q, V) using Mayer-Vietoris for Q = C_8 * C_8.
    dim_H0_Q_V = dim_H0_C8_V
    dim_H1_Q_V = dim_V - dim_H0_Q_V
    
    print("Step 3: Compute dimensions of H^p(Q, V), where Q=C_8*C_8.")
    print(f"dim H^0(Q, V) = {dim_H0_Q_V}")
    print(f"dim H^1(Q, V) = dim(V) - dim H^0(Q,V) = {dim_V} - {dim_H0_Q_V} = {dim_H1_Q_V}")
    print(f"dim H^p(Q, V) = 0 for p >= 2.")
    print("-" * 30)
    
    # Step 4: Use LHS spectral sequence to find dim H^2(G,M).
    # E_2^{p,q} = H^p(Q, H^q(Z,M)). We need sum of dim E_inf^{p,q} for p+q=2.
    # Relevant E_2 terms are E_2^{2,0}, E_2^{1,1}, E_2^{0,2}.
    # The spectral sequence collapses at E_2 for this diagonal, so dim H^2(G,M) = dim E_2^{1,1}.
    dim_E2_2_0 = 0 # This is H^2(Q,V)
    dim_E2_0_2 = 0 # This is H^0(Q,0)
    dim_E2_1_1 = dim_H1_Q_V
    dim_H2_G_M = dim_E2_2_0 + dim_E2_1_1 + dim_E2_0_2

    print("Step 4: Combine results using the LHS spectral sequence.")
    print("The E_2 page has terms E_2^{p,q} = H^p(Q, H^q(Z,M)).")
    print("For p+q=2, the only non-zero term is E_2^{1,1} = H^1(Q, H^1(Z,M)).")
    print(f"dim E_2^{2,0} = dim H^2(Q, H^0(Z,M)) = {dim_E2_2_0}")
    print(f"dim E_2^{1,1} = dim H^1(Q, H^1(Z,M)) = {dim_E2_1_1}")
    print(f"dim E_2^{0,2} = dim H^0(Q, H^2(Z,M)) = {dim_E2_0_2}")
    
    print("\nFinal calculation:")
    print(f"dim H^2(G,M) = dim E_infinity^{2,0} + dim E_infinity^{1,1} + dim E_infinity^{0,2}")
    print(f"               = {dim_E2_2_0} + {dim_E2_1_1} + {dim_E2_0_2} = {dim_H2_G_M}")

solve()