import math

def solve():
    """
    Solves the problem by calculating f(alpha_p, beta_p, gamma_p) mod p.
    """

    # --- Step 1: Define the function f(a,b,c) and helpers ---
    # f(a,b,c) is the number of paths from (0,0,0) to (a,b,c) using steps
    # (1,0,0), (0,2,0), and (0,0,3). This can be calculated with a
    # multinomial coefficient.
    
    memo_factorial = {}
    def factorial(n):
        if n in memo_factorial:
            return memo_factorial[n]
        if n == 0:
            return 1
        res = 1
        for i in range(1, n + 1):
            res *= i
        memo_factorial[n] = res
        return res

    def f(a, b, c):
        # f(a,b,c) = 0 if b is not even or c is not a multiple of 3.
        if b % 2 != 0 or c % 3 != 0:
            return 0
        
        n1 = a
        n2 = b // 2
        n3 = c // 3
        
        # Multinomial coefficient: (n1+n2+n3)! / (n1! * n2! * n3!)
        numerator = factorial(n1 + n2 + n3)
        denominator = factorial(n1) * factorial(n2) * factorial(n3)
        return numerator // denominator

    # --- Step 2: Use Lucas's Theorem Analogue ---
    # The problem reduces to computing a product of f(a,b,c) for the base-p digits.
    # The coefficient vectors (A_k, B_k, C_k) from the base-p expansion of 
    # (alpha_p, beta_p, gamma_p) repeat in a cycle of 3.
    coeffs = [
        (1, 8, 3),   # for k = 3i
        (3, 4, 9),   # for k = 3i+1
        (4, 4, 12)   # for k = 3i+2
    ]

    # --- Step 3: Calculate the component f values ---
    f_values = [f(a, b, c) for a, b, c in coeffs]
    f0, f1, f2 = f_values

    print("Step-by-step Calculation Plan:")
    print("1. The problem simplifies using a property similar to Lucas's Theorem.")
    print("   f(alpha_p, beta_p, gamma_p) mod p = Product[f(A_k, B_k, C_k)] mod p")
    print("   where (A_k, B_k, C_k) are the base-p digits of the inputs.")
    print("\n2. The base-p digits repeat in a cycle of 3. We calculate f for each unique vector:")
    print(f"   f{coeffs[0]} = {f0}")
    print(f"   f{coeffs[1]} = {f1}")
    print(f"   f{coeffs[2]} = {f2}")

    # --- Step 4: Calculate the product P ---
    P = f0 * f1 * f2
    print(f"\n3. Let P be the product of these values: P = {f0} * {f1} * {f2} = {P}")

    # --- Step 5 & 6: Simplify exponent and compute Legendre Symbol ---
    # The final value is P^((3p+1)/2) mod p.
    # This simplifies to P^2 * (P/p) mod p, where (P/p) is the Legendre symbol.
    # (P/p) = (52920000/p) = ( (2^6 * 3^3 * 5^4 * 7^2) / p ) = (3/p).
    # Using quadratic reciprocity for p = 2^127-1:
    # p mod 3 = 1, p mod 4 = 3.
    # (3/p) = (p/3) * (-1)^((p-1)/2 * (3-1)/2) = 1 * (-1)^((p-1)/2) = -1.
    legendre_symbol = -1
    print("\n4. The expression simplifies to P^2 * (P/p) mod p.")
    print(f"   The Legendre symbol (P/p) is calculated to be {legendre_symbol}.")

    # --- Step 7: Final Calculation ---
    # Result = P^2 * (-1) mod p = -P^2 mod p.
    P_squared = P**2
    p = 2**127 - 1
    final_result = p - P_squared

    print(f"\n5. The final result is -P^2 mod p.")
    print("\nFinal Equation:")
    print(f"f(alpha_p, beta_p, gamma_p) mod p = (2^127 - 1) - {P}^2")
    print(f"f(alpha_p, beta_p, gamma_p) mod p = {p} - {P_squared}")
    print(f"f(alpha_p, beta_p, gamma_p) mod p = {final_result}")

solve()
<<<170141183460469231731687303715881306199>>>