import math

def gaussian_binomial_coefficient(n, k, q):
    """Calculates the Gaussian binomial coefficient."""
    if k < 0 or k > n:
        return 0
    num = 1
    for i in range(n - k + 1, n + 1):
        num *= (q**i - 1)
    den = 1
    for i in range(1, k + 1):
        den *= (q**i - 1)
    return num // den

def proportion_calculation():
    """
    Calculates the proportion of irreducible (3,2)-stingray duos.
    """
    q = 4  # size of the finite field
    d = 5  # dimension of the vector space
    e1 = 3 # dimension of U1
    f1 = d - e1 # dimension of F1
    e2 = 2 # dimension of U2
    f2 = d - e2 # dimension of F2
    
    # We want to calculate P(irreducible | duo).
    # This is 1 - P(reducible | duo).
    # A duo is reducible if F1^F2 != {0} (E1), U1=F2 (E2), or U2=F1 (E3).
    # These events are mutually exclusive for a duo.
    # P(red|duo) = P(E1|C) + P(E2|C) + P(E3|C), where C is U1^U2 = {0}.
    
    print("Step 1: Calculate probability of the duo condition C: U1 intersect U2 = {0}")
    # This is the probability that a random 2-space (U2) is disjoint from a fixed 3-space (U1).
    # Number of 2-spaces disjoint from a fixed 3-space U1:
    # N_disjoint = q^(e1*e2) * GBC(d-e1, e2, q)
    N_U2_disjoint_U1 = q**(e1 * e2) * gaussian_binomial_coefficient(d - e1, e2, q)
    # Total number of 2-spaces:
    N_U2_total = gaussian_binomial_coefficient(d, e2, q)
    P_C = N_U2_disjoint_U1 / N_U2_total
    print(f"Total number of 2-dim subspaces in F_4^5: {N_U2_total}")
    print(f"Number of 2-dim subspaces disjoint from a fixed 3-dim subspace: {N_U2_disjoint_U1}")
    print(f"P(C) = P(U1 intersect U2 = {{0}}) = {N_U2_disjoint_U1}/{N_U2_total}\n")

    print("Step 2: Calculate probability of reducibility from condition (1): F1 intersect F2 != {0}")
    # P(E1|C) = P(F1^F2 != {0} | U1^U2 = {0}).
    # We assume independence of these subspace events. So P(E1|C) approx P(E1).
    # P(E1) is probability that a random 3-space (F2) intersects a fixed 2-space (F1).
    # Total number of 3-spaces:
    N_F2_total = gaussian_binomial_coefficient(d, f2, q)
    # Number of 3-spaces disjoint from a fixed 2-space F1:
    N_F2_disjoint_F1 = q**(f1 * f2) * gaussian_binomial_coefficient(d - f1, f2, q)
    # P(F1 intersect F2 = {0})
    P_F1_F2_disjoint = N_F2_disjoint_F1 / N_F2_total
    P_E1 = 1 - P_F1_F2_disjoint
    print(f"Total number of 3-dim subspaces: {N_F2_total}")
    print(f"Number of 3-dim subspaces disjoint from a fixed 2-dim subspace: {N_F2_disjoint_F1}")
    print(f"P(E1) = P(F1 intersect F2 != {{0}}) = 1 - {N_F2_disjoint_F1}/{N_F2_total} = {(N_F2_total - N_F2_disjoint_F1)}/{N_F2_total}\n")
    P_E1_given_C = P_E1 # by independence assumption
    
    print("Step 3: Calculate probability of reducibility from condition (2): U1 = F2")
    # P(E2|C) = P(U1=F2 | C) = P(U1=F2 and C) / P(C).
    # If U1=F2, then V = U2+F2 = U2+U1, which implies U1 intersect U2 = {0}.
    # So P(U1=F2 and C) = P(U1=F2).
    # P(U1=F2) is probability a random 3-space is a specific 3-space.
    P_E2 = 1 / N_F2_total
    P_E2_given_C = P_E2 / P_C
    print(f"P(E2) = P(U1 = F2) = 1/{N_F2_total}")
    print(f"P(E2|C) = P(E2)/P(C) = (1/{N_F2_total}) / ({N_U2_disjoint_U1}/{N_U2_total}) = {N_U2_total}/({N_F2_total}*{N_U2_disjoint_U1}) = 1/{N_U2_disjoint_U1}\n")

    print("Step 4: Calculate probability of reducibility from condition (3): U2 = F1")
    # P(E3|C) = P(U2=F1 | C).
    # If U2=F1, then V = U1+F1 = U1+U2, which implies U1 intersect U2 = {0}.
    # So P(U2=F1 and C) = P(U2=F1).
    # P(U2=F1) is probability a random 2-space is a specific 2-space.
    P_E3 = 1 / N_U2_total
    P_E3_given_C = P_E3 / P_C
    print(f"P(E3) = P(U2 = F1) = 1/{N_U2_total}")
    # Note: N_U2_total == N_F2_total because GBC(n,k)=GBC(n,n-k)
    print(f"P(E3|C) = P(E3)/P(C) = (1/{N_U2_total}) / ({N_U2_disjoint_U1}/{N_U2_total}) = 1/{N_U2_disjoint_U1}\n")
    
    print("Step 5: Calculate total probability of being reducible and irreducible.")
    # The three conditions for reducibility are mutually exclusive.
    P_red_given_C_num = (N_F2_total - N_F2_disjoint_F1) * N_U2_disjoint_U1 + N_F2_total + N_U2_total
    P_red_given_C_den = N_F2_total * N_U2_disjoint_U1
    
    p_e1_c = (N_F2_total - N_F2_disjoint_F1) / N_F2_total
    p_e2_c = 1 / N_U2_disjoint_U1
    p_e3_c = 1 / N_U2_disjoint_U1
    
    P_red_given_C = p_e1_c + p_e2_c + p_e3_c
    
    # Using fractions for precision
    from fractions import Fraction
    frac_p_e1_c = Fraction(N_F2_total - N_F2_disjoint_F1, N_F2_total)
    frac_p_e2_c = Fraction(1, N_U2_disjoint_U1)
    frac_p_e3_c = Fraction(1, N_U2_disjoint_U1)

    frac_P_red = frac_p_e1_c + frac_p_e2_c + frac_p_e3_c
    frac_P_irr = 1 - frac_P_red

    print(f"P(reducible|duo) = P(E1|C) + P(E2|C) + P(E3|C)")
    print(f"= {N_F2_total - N_F2_disjoint_F1}/{N_F2_total} + 1/{N_U2_disjoint_U1} + 1/{N_U2_disjoint_U1}")
    print(f"= {frac_p_e1_c} + {frac_p_e2_c} + {frac_p_e3_c} = {frac_P_red}")

    print(f"\nP(irreducible|duo) = 1 - P(reducible|duo) = 1 - {frac_P_red} = {frac_P_irr}")

    print("\nFinal Answer for (c):")
    # We output the simplified fraction and its decimal value
    print(f"The proportion is {frac_P_irr.numerator}/{frac_P_irr.denominator}, which is approximately {float(frac_P_irr)}.")
    
proportion_calculation()