import math

def get_divisors(n):
    """Returns the set of divisors of n."""
    if n == 0:
        return []
    divs = {1}
    if n > 1:
        divs.add(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return list(divs)

def sigma(k, n):
    """Computes sigma_k(n), the sum of the k-th powers of the divisors of n."""
    if n <= 0:
        return 0
    if n == 1:
        return 1
    divisors = get_divisors(n)
    return sum(d**k for d in divisors)

def get_e4_coeffs(num_coeffs):
    """Generates coefficients of the E_4(z) q-expansion."""
    coeffs = [1]
    for n in range(1, num_coeffs + 1):
        coeffs.append(240 * sigma(3, n))
    return coeffs

def get_e8_coeffs(num_coeffs):
    """Generates coefficients of the E_8(z) q-expansion."""
    coeffs = [1]
    for n in range(1, num_coeffs + 1):
        coeffs.append(480 * sigma(7, n))
    return coeffs

def multiply_series(A, B, max_n):
    """Multiplies two q-series A and B up to the q^max_n term."""
    C = [0] * (max_n + 1)
    for i in range(max_n + 1):
        for j in range(i + 1):
            if j < len(A) and (i - j) < len(B):
                C[i] += A[j] * B[i - j]
    return C

def main():
    """
    Main function to calculate the sum of the first three non-zero coefficients
    of the normalized cusp form.
    """
    # The problem asks for a cusp form in the subspace spanned by g1, g2, g3.
    # g1 = E4^2, g2 = E4*E4(2z), g3 = E4(2z)^2
    # Theory shows the unique cusp form is proportional to g1 - g2.

    # We need coefficients for g1 and g2.
    # g1 = E4^2 = E8.
    num_coeffs = 10
    g1_coeffs = get_e8_coeffs(num_coeffs)

    # For g2 = E4(z) * E4(2z)
    e4_coeffs = get_e4_coeffs(num_coeffs)
    e4_2z_coeffs = [0] * (num_coeffs + 1)
    for i in range(len(e4_coeffs)):
        if 2 * i < len(e4_2z_coeffs):
            e4_2z_coeffs[2 * i] = e4_coeffs[i]

    g2_coeffs = multiply_series(e4_coeffs, e4_2z_coeffs, num_coeffs)

    # The unnormalized cusp form is f_un = g1 - g2.
    f_un_coeffs = [g1_coeffs[i] - g2_coeffs[i] for i in range(num_coeffs + 1)]

    # Find the first non-zero coefficient to normalize.
    # This will be the coefficient of q^1.
    norm_factor = f_un_coeffs[1]
    
    # The normalized coefficients are a_n = f_un_coeffs[n] / norm_factor
    # Find the first three non-zero coefficients of the normalized form.
    
    a_coeffs_normalized = []
    for i in range(1, num_coeffs + 1):
        coeff = f_un_coeffs[i] / norm_factor
        # Since coeffs will be integers, round to avoid float precision issues.
        a_coeffs_normalized.append(round(coeff))

    # The problem implies coefficients are non-zero starting from a_1.
    first_three_coeffs = a_coeffs_normalized[:3]
    total_sum = sum(first_three_coeffs)
    
    c1, c2, c3 = first_three_coeffs
    
    print(f"The first three non-zero coefficients of the normalized cusp form f are: {c1}, {c2}, and {c3}.")
    print(f"The sum is: {c1} + {c2} + {c3} = {total_sum}")


main()