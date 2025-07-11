import math

def solve_m():
    """
    Solves for the smallest positive integer M based on the problem description.
    """
    N = 100000

    # The problem is equivalent to finding the smallest S such that for any a_i in [0,1],
    # there is a choice of signs sigma_i such that the range of partial sums P_k = sum_{i=1 to k} sigma_i a_i
    # satisfies max(P_k) - min(P_k) <= 2S.
    # Let C_N be the supremum of this minimal range over all sequences a_i.
    # The condition is 2S >= C_N.
    # For large even N, it is a known (though advanced) result from discrepancy theory that
    # C_N = 2 - 2/N.
    
    # Step 1: Calculate C_N for N = 100000
    C_N = 2 - 2 / N
    
    print(f"Let N = {N}.")
    print(f"The maximum minimal range of the partial sums, C_N, is given by the formula C_N = 2 - 2/N.")
    print(f"C_N = 2 - 2/{N} = {C_N}")
    
    # Step 2: Use the condition on S = 1 - 10^(-M)
    # 2*S >= C_N
    # 2 * (1 - 10**(-M)) >= C_N
    S_required = C_N / 2
    
    print(f"\nThe condition on S is 2*S >= C_N, which means S >= C_N / 2.")
    print(f"S >= {C_N} / 2 = {S_required}")
    
    # Step 3: Solve the inequality for M
    # 1 - 10**(-M) >= S_required
    # 1 - S_required >= 10**(-M)
    # log10(1 - S_required) >= -M
    # M >= -log10(1 - S_required)
    
    # S_required = (2 - 2/N) / 2 = 1 - 1/N
    # 1 - S_required = 1 - (1 - 1/N) = 1/N
    one_minus_S_required = 1 / N
    
    print(f"\nWe have S = 1 - 10**(-M). So we must solve the inequality:")
    print(f"1 - 10**(-M) >= {S_required}")
    print(f"1 - {S_required} >= 10**(-M)")
    print(f"{1 - S_required:.6f} >= 10**(-M)")
    
    # We can write 1/N in scientific notation
    print(f"Since 1/N = 1/{N} = {one_minus_S_required}, the inequality is:")
    print(f"10**(-5) >= 10**(-M)")
    
    # Taking log base 10 of both sides
    # -5 >= -M
    # M >= 5
    
    print(f"\nTaking the base-10 logarithm of both sides gives:")
    print(f"log10(10**(-5)) >= log10(10**(-M))")
    print(f"-5 >= -M")
    print(f"Multiplying by -1 and flipping the inequality sign gives:")
    print(f"5 <= M")
    
    # The smallest positive integer M is ceil(-log10(1/N))
    M = math.ceil(-math.log10(one_minus_S_required))
    
    print(f"\nThe smallest positive integer M that satisfies M >= 5 is {M}.")

solve_m()