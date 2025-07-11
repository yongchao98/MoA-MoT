import math

def sigma(k, n):
    """
    Calculates the sum of the k-th powers of the divisors of n.
    sigma_k(n) = sum_{d|n} d^k
    """
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

def get_e4_coeffs(num_coeffs):
    """
    Computes the first `num_coeffs` coefficients of the q-expansion of E_4(z).
    E_4(z) = 1 + 240 * sum_{n=1 to inf} sigma_3(n) * q^n
    """
    coeffs = [1]
    # The constant factor for E4 is -2*4/B_4 = -8/(-1/30) = 240
    C4 = 240
    for n in range(1, num_coeffs):
        coeffs.append(C4 * sigma(3, n))
    return coeffs

def poly_multiply(p1, p2, max_deg):
    """
    Multiplies two polynomials represented as lists of coefficients.
    """
    n1 = len(p1)
    n2 = len(p2)
    result = [0] * max_deg
    for i in range(min(n1, max_deg)):
        for j in range(min(n2, max_deg - i)):
            result[i+j] += p1[i] * p2[j]
    return result

def main():
    """
    Main function to calculate the sum of the coefficients.
    """
    # We need at least the first 4 coefficients (for q^0 to q^3)
    NUM_COEFFS = 5

    # 1. Get coefficients for E_4(z)
    e4_coeffs = get_e4_coeffs(NUM_COEFFS)
    # e4_coeffs: [a_0, a_1, a_2, a_3, a_4, ...]

    # 2. Get coefficients for F(z) = E_4(2z)
    # The q-expansion is E_4(q^2), so coeffs are spaced out.
    f_coeffs = [0] * NUM_COEFFS
    for i in range(NUM_COEFFS):
        if i % 2 == 0 and i // 2 < len(e4_coeffs):
            f_coeffs[i] = e4_coeffs[i//2]

    # 3. Get coefficients for the basis g1=E_4^2 and g3=F^2
    g1_coeffs = poly_multiply(e4_coeffs, e4_coeffs, NUM_COEFFS)
    g3_coeffs = poly_multiply(f_coeffs, f_coeffs, NUM_COEFFS)

    # 4. Construct the un-normalized cusp form f_un = g1 - g3
    f_un_coeffs = [g1_coeffs[i] - g3_coeffs[i] for i in range(NUM_COEFFS)]

    # 5. Normalize the cusp form f by dividing by the first coefficient a_1
    # The constant term f_un_coeffs[0] is 0, as expected.
    # The first coefficient is f_un_coeffs[1].
    first_coeff = f_un_coeffs[1]
    if first_coeff == 0:
        print("Error: The first coefficient is zero, cannot normalize.")
        return

    normalized_coeffs = [c / first_coeff for c in f_un_coeffs]

    # 6. Find the sum of the first three non-zero coefficients.
    # These are the coefficients of q^1, q^2, and q^3.
    c1 = int(round(normalized_coeffs[1]))
    c2 = int(round(normalized_coeffs[2]))
    c3 = int(round(normalized_coeffs[3]))

    total_sum = c1 + c2 + c3
    
    # 7. Print the final equation and sum
    print("The q-expansion of the normalized cusp form f(z) begins:")
    print(f"f(q) = {c1}q + {c2}q^2 + {c3}q^3 + ...")
    print("\nThe sum of the first three non-zero coefficients is:")
    print(f"{c1} + {c2} + {c3} = {total_sum}")


if __name__ == "__main__":
    main()
<<<2317>>>