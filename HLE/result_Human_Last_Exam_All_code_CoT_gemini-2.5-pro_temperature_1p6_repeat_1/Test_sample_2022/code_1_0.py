import math

def solve():
    """
    Calculates the value of F(n) for two prime numbers.
    """

    # The two prime numbers
    p1 = 80039
    p2 = 80077

    # From the problem analysis, n simplifies to 6 mod p.
    # We need to calculate S(6).
    # S(6) = sum_{k=0 to 3} C(6-2k, k) * (-1/8)^k
    # S(6) = C(6,0)*(-1/8)^0 + C(4,1)*(-1/8)^1 + C(2,2)*(-1/8)^2 + C(0,3)*(-1/8)^3
    
    # Calculate terms of S(6)
    term0_num = math.comb(6, 0)
    term0_den = 1
    
    term1_num = math.comb(4, 1) * (-1)
    term1_den = 8

    term2_num = math.comb(2, 2)
    term2_den = 64

    # S(6) = 1 - 1/2 + 1/64 = 33/64
    a = 33
    b = 64
    
    # Print the equation as requested
    print("The problem requires calculating F(n) mod p, where n = p^5+2p^4-19p^3-3p^2+16p+6.")
    print("First, we simplify n mod p, which results in n = 6.")
    print("So we need to calculate F(6) = S(6) mod p.")
    print("\nCalculating S(6):")
    print(f"S(6) = (C(6,0) * (-1/8)^0) + (C(4,1) * (-1/8)^1) + (C(2,2) * (-1/8)^2)")
    print(f"S(6) = (1 * 1) + (4 * -1/8) + (1 * 1/64)")
    print(f"S(6) = 1 - 1/2 + 1/64")
    print(f"S(6) = {a}/{b}")
    
    # Calculate F(6) for p1
    # F(6) mod p1 = (a * b^-1) mod p1
    b_inv_p1 = pow(b, -1, p1)
    result1 = (a * b_inv_p1) % p1
    
    # Calculate F(6) for p2
    # F(6) mod p2 = (a * b^-1) mod p2
    b_inv_p2 = pow(b, -1, p2)
    result2 = (a * b_inv_p2) % p2
    
    print(f"\nFor p = {p1}, the value is: {result1}")
    print(f"For p = {p2}, the value is: {result2}")
    
    # Final answer format
    print("\nThe final answers separated by a comma are:")
    print(f"{result1},{result2}")

solve()