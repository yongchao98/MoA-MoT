import math

def solve():
    """
    Solves for the smallest positive integer M based on the problem description.
    """
    
    # Step 1: Find the largest integer k such that the harmonic number H_k is at most 2.
    # This value, H_k, corresponds to L.
    k = 0
    harmonic_sum = 0.0
    k_val = 0
    
    while True:
        k_val += 1
        next_sum = harmonic_sum + 1.0 / k_val
        if next_sum > 2:
            break
        harmonic_sum = next_sum
        k = k_val

    L = harmonic_sum
    
    # To display L as a fraction, we compute the numerator and denominator for H_k.
    # For k=3, H_3 = 1/1 + 1/2 + 1/3 = 6/6 + 3/6 + 2/6 = 11/6.
    if k == 3:
      L_num = 11
      L_den = 6
    else: # A general calculation for H_k as a fraction
      den = 1
      for i in range(1, k + 1):
        den = den * i // math.gcd(den, i)
      num = 0
      for i in range(1, k + 1):
        num += den // i
      common = math.gcd(num, den)
      L_num = num // common
      L_den = den // common

    # Step 2: The condition is L <= 2 * (1 - 10**-M).
    # We solve for M.
    # L/2 <= 1 - 10**-M
    # 10**-M <= 1 - L/2
    # -M <= log10(1 - L/2)
    # M >= -log10(1 - L/2)
    # M >= log10(2 / (2 - L))
    
    m_threshold = math.log10(2 / (2 - L))
    
    # M must be the smallest positive integer greater than or equal to this threshold.
    M = math.ceil(m_threshold)

    print("The problem requires finding the value L, which is the maximum possible minimal range of the sums.")
    print(f"This value L is the largest harmonic number H_k <= 2, which is H_{k} = {L_num}/{L_den}.")
    print("\nThe inequality to be solved is:")
    print(f"{L_num}/{L_den} <= 2 * (1 - 10**-{M})")
    
    print("\nSolving for M:")
    print(f"M >= log10(2 / (2 - {L_num}/{L_den}))")
    print(f"M >= log10(12)")
    print(f"M >= {m_threshold:.4f}")
    
    print(f"\nThe smallest positive integer M is {M}.")
    
solve()