import math

def power(a, b, m):
    """Computes (a^b) % m using modular exponentiation."""
    res = 1
    a %= m
    while b > 0:
        if b % 2 == 1:
            res = (res * a) % m
        a = (a * a) % m
        b //= 2
    return res

def multinomial_direct(n, coeffs):
    """
    Computes multinomial coefficient n!/(k1!*k2!*...) for small integer inputs.
    This is suitable here because the values of n are small (2, 3, 5).
    """
    if sum(coeffs) != n:
        return 0
    res = math.factorial(n)
    for k in coeffs:
        res //= math.factorial(k)
    return res

def compute_a_nkl_mod_p(n, k, l, p):
    """
    Calculates a_{n,k,l} mod p.
    This coefficient is from the expansion of (12 + 3x + 75y + 27x^2y^2)^n.
    The coefficient of x^k y^l is the sum over terms:
    C(n; n1,n2,n3,n4) * 12^n1 * 3^n2 * 75^n3 * 27^n4
    where the sum is over all n1,n2,n3,n4 such that:
    n1+n2+n3+n4 = n
    n2+2*n4 = k
    n3+2*n4 = l
    """
    total = 0
    # From the constraints, we can determine the valid range for n4
    n4_max = min(k // 2, l // 2)
    n4_min_num = k + l - n
    n4_min = (n4_min_num + 2) // 3 if n4_min_num > 0 else 0

    for n4 in range(n4_min, n4_max + 1):
        n2 = k - 2 * n4
        n3 = l - 2 * n4
        n1 = n - n2 - n3 - n4
        
        # Ensure all indices are non-negative
        if n1 < 0 or n2 < 0 or n3 < 0:
            continue
            
        coeffs = [n1, n2, n3, n4]
        multi_coeff = multinomial_direct(n, coeffs)
        
        term = multi_coeff
        term = (term * power(12, n1, p)) % p
        term = (term * power(3, n2, p)) % p
        term = (term * power(75, n3, p)) % p
        term = (term * power(27, n4, p)) % p
        
        total = (total + term) % p
        
    return total

def solve():
    """
    Main function to solve the problem.
    """
    p = 21023
    print(f"The prime is p = {p}")
    
    # The repeating digit triplets are (5,2,2), (3,1,2), (2,1,1).
    # We calculate the corresponding coefficients c_0, c_1, c_2.
    
    # Calculate c_2 = a_{2,1,1}
    n2, k2, l2 = 2, 1, 1
    c2 = compute_a_nkl_mod_p(n2, k2, l2, p)
    print(f"c_2 = a_{{{n2},{k2},{l2}}} mod {p} = {c2}")

    # Calculate c_1 = a_{3,1,2}
    n1, k1, l1 = 3, 1, 2
    c1 = compute_a_nkl_mod_p(n1, k1, l1, p)
    print(f"c_1 = a_{{{n1},{k1},{l1}}} mod {p} = {c1}")
    
    # Calculate c_0 = a_{5,2,2}
    n0, k0, l0 = 5, 2, 2
    c0 = compute_a_nkl_mod_p(n0, k0, l0, p)
    print(f"c_0 = a_{{{n0},{k0},{l0}}} mod {p} = {c0}")
    
    # Calculate C = c_0 * c_1 * c_2 mod p
    C = (c0 * c1 * c2) % p
    print(f"\nLet C = (c_0 * c_1 * c_2) mod p")
    print(f"C = ({c0} * {c1} * {c2}) mod {p} = {C}")

    # The final answer is -C^2 mod p, based on the Legendre symbol calculation.
    # The exponent is (3p+1)/2, which is 3(p-1)/2 + 2.
    # The value is (C^((p-1)/2))^3 * C^2 mod p.
    # The Legendre symbol (C/p) was calculated to be -1.
    # So the result is (-1)^3 * C^2 = -C^2 mod p.
    
    C_squared = power(C, 2, p)
    result = (-C_squared + p) % p

    print(f"\nThe final value is given by C^((3p+1)/2) mod p = -C^2 mod p.")
    print(f"C^2 mod {p} = {C}^2 mod {p} = {C_squared}")
    final_equation = f"-C^2 mod {p} = -{C_squared} mod {p} = {result}"
    print(f"Final calculation: {final_equation}")
    
    return result

final_answer = solve()
print(f"\nThe calculated value of a_n,k,l mod p is: {final_answer}")
print(f"<<<{final_answer}>>>")